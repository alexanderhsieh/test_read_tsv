version 1.0


###########################################################################
## WORKFLOW DEFINITION
###########################################################################
workflow process_sample_table {
	input {
		File table
		
		File ref_fasta
		File ref_fasta_index
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
		Array[String] selected_readgroup_id_columns = select_all(rg_id_columns)

		call merge_trio_gvcf {
			input:
				sample_id = sample_id,
				trio_gvcf_array = selected_gvcf_columns,
				trio_gvcf_index_array = selected_gvcf_index_columns,
				trio_readgroup_ids = selected_readgroup_id_columns,
				ref_fasta = ref_fasta,
				ref_fasta_index = ref_fasta_index
		}

	}



	output {
		Array[File] merged_trio_gvcfs = merge_trio_gvcf.out_gvcf
		Array[File] merged_trio_gvcf_index = merge_trio_gvcf.out_gvcf_index
	}



}



###########################################################################
## TASK DEFINITIONS
###########################################################################
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

task merge_trio_gvcf {
	input{
		String sample_id

		Array[File] trio_gvcf_array
		Array[File] trio_gvcf_index_array
		Array[String] trio_readgroup_ids

		File ref_fasta
		File ref_fasta_index
	}


	String outfname = "~{sample_id}.TRIO.g.vcf.gz"

	command {

		bcftools merge -g ~{ref_fasta} -l ~{write_lines(trio_gvcf_array)} -o ~{outfname} -O z -m all

		tabix -p vcf ~{outfname}

		cat ~{write_lines(trio_readgroup_ids)} > id_list.txt
		sed -n '1p' id_list.txt = pb_id.txt
		sed -n '2p' id_list.txt = fa_id.txt
		sed -n '3p' id_list.txt = mo_id.txt


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

		File out_pb_id = "pb_id.txt"
		File out_fa_id = "fa_id.txt"
		File out_mo_id = "mo_id.txt"
	}

}

