version 1.0


###########################################################################
## WORKFLOW DEFINITION
###########################################################################
workflow gvcf_to_denovo {
	input {
		File table
		
		File ref_fasta
		File ref_fasta_index

		File dn_script 

		Float pb_min_vaf
		Int par_min_dp
		Int par_max_alt

		String output_prefix
		String output_suffix
	}

	call read_table {
		input:
			table = table
	}

	scatter (i in range(length(read_table.out))) {

		Int n_cols = length(read_table.out[0])

		## parse columns containing gvcf google bucket paths 
		scatter (j in range(n_cols)) {
			if (j >=1 && j<=3) {
				File gvcf_columns = read_table.out[i][j]
			}
		}

		## parse columns containing gvcf index google bucket paths
		scatter (j in range(n_cols)) {
			if (j >=4 && j<=6) {
				File gvcf_index_columns = read_table.out[i][j]
			}
		}

		## parse columns containing gvcf readgroup ids
		scatter (j in range(n_cols)) {
			if (j >=7 && j<=9) {
				File rg_id_columns = read_table.out[i][j]
			}
		}

		String sample_id = read_table.out[i][0]
		Array[File] selected_gvcf_columns = select_all(gvcf_columns)
		Array[File] selected_gvcf_index_columns = select_all(gvcf_index_columns)
		#Array[String] selected_readgroup_id_columns = select_all(rg_id_columns)

		call merge_trio_gvcf {
			input:
				sample_id = sample_id,
				trio_gvcf_array = selected_gvcf_columns,
				trio_gvcf_index_array = selected_gvcf_index_columns,
				ref_fasta = ref_fasta,
				ref_fasta_index = ref_fasta_index
		}

		call call_denovos {
			input:
				script = dn_script,
				sample_id = sample_id,
				trio_readgroup_ids = merge_trio_gvcf.rg_ids,
				gvcf = merge_trio_gvcf.out_gvcf,
				pb_min_vaf = pb_min_vaf,
				par_max_alt = par_max_alt,
				par_min_dp = par_min_dp,
				output_suffix = output_suffix
		}

	}

	call gather_shards {
		input:
			shards = call_denovos.out,
			headers = call_denovos.head,
			output_prefix = output_prefix,
			output_suffix = output_suffix
	}


	output {
		File raw_denovos = gather_shards.out
	}



}



###########################################################################
## TASK DEFINITIONS
###########################################################################
## Reads in sample_attributes.no_header.tsv to enable coercion from 
## google bucket_path (String) to corresponding gvcf file (File)
## Note: requires specific 10-column format:
##		 sample_id, sample_gvcf, father_gvcf, mother_gvcf, 
##		 sample_gvcf_index, father_gvcf_index, mother_gvcf_index, 
##		 sample_rg_id, father_rg_id, mother_rg_id
task read_table {
	input{
		File table
	}

	command { 
		echo "reading table" 
	}

	runtime{
		docker: "ubuntu:latest"
	}

	output{
		Array[Array[String]] out = read_tsv(table)
	}

}


## Runs bcftools merge on the trio gvcf array supplied in read_table()
## Note: requires array of gvcf indices to be present as well
task merge_trio_gvcf {
	input{
		String sample_id

		Array[File] trio_gvcf_array
		Array[File] trio_gvcf_index_array

		File ref_fasta
		File ref_fasta_index
	}

	String outfname = "~{sample_id}.TRIO.g.vcf.gz"

	command {

		bcftools merge -g ~{ref_fasta} -l ~{write_lines(trio_gvcf_array)} -o ~{outfname} -O z -m all

		tabix -p vcf ~{outfname}

		bcftools query -l ~{outfname} > "ids.txt"


	}

	runtime {
		docker: "gatksv/sv-base-mini:cbb1fc"
		memory: "8G"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File out_gvcf = "~{outfname}"
		File out_gvcf_index = "~{outfname}.tbi"
		Array[String] rg_ids = read_lines("ids.txt")
	}

}

## Calls denovos from bcftools merged gvcf containing proband, father, mother
## Note: readgroup IDs in the GVCF may differ from that in the PED or Sample Map
##		 here, we parse from the GVCF filename upstream when creating sample_attributes.tsv 
task call_denovos {
	input {
		File script

		String sample_id

		Array[String] trio_readgroup_ids

		File gvcf

		Float pb_min_vaf
		Int par_max_alt
		Int par_min_dp

		String output_suffix
	}

	String pb_id = trio_readgroup_ids[0]
	String fa_id = trio_readgroup_ids[1]
	String mo_id = trio_readgroup_ids[2]

	String output_file = "~{sample_id}~{output_suffix}"

	command {

		echo "~{pb_id} ~{fa_id} ~{mo_id}"

		python3 ${script} -s ~{pb_id} -f ~{fa_id} -m ~{mo_id} -g ~{gvcf} -x ~{pb_min_vaf} -y ~{par_max_alt} -z ~{par_min_dp} -o ~{output_file}

		grep "^id" ~{output_file} > "header.txt"
	}

	runtime {
		docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File out = "~{output_file}"
		File head = "header.txt"
	}
}

#Gathers shards of de novo calls into a single callset
task gather_shards {
	input {
		Array[File] shards 
		Array[File] headers
		String output_prefix
		String output_suffix
	}

	File header = headers[0]
	String output_file = "~{output_prefix}~{output_suffix}"

	command {

		while read file; do
		cat $file | grep -v "^id" >> "tmp.cat.txt"
		done < ~{write_lines(shards)};


		## CHECK THAT DE NOVO CALLSETS ARE NOT EMPTY
		N_LINES=`wc -l "tmp.cat.txt"`
		echo "tmp cat: $N_LINES lines"

		if [ $N_LINES -le 1 ]; then
			echo "EMPTY OUTPUT DE NOVO CALLSET"; exit $ERRCODE; 
		else
			(cat "~{header}" "tmp.cat.txt") > "~{output_file}"
		fi

		



	}

	runtime {
		docker: "ubuntu:latest"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File out = "~{output_file}"
	}
}


