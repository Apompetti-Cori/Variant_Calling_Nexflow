nextflow -log $(pwd)/logs/nextflow.log \
run $(find $(pwd) -maxdepth 2 -type f -name "main.nf") \
-resume \
-entry somatic_variants \
-with-report $(pwd)/logs/run-report-$(date +"%Y-%m-%d_%I:%M:%p").html \
-with-trace $(pwd)/logs/run-trace-$(date +"%Y-%m-%d_%I:%M:%p").txt \
-with-dag $(pwd)/logs/run-flowchart.mmd \
--genome "hg38" --multiqc_report_title "Test Title" --sample_table $(find $(pwd) -type f -name "sample_table.csv")