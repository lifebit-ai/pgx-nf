# pgx-nf

## Quick Start
The tool(s) can be run on:
* [the command line](#running-on-the-command-line)
* [Deploit](#deploit) (recommended)

## Running on the command line

### Dependencies 
[Nextflow](https://www.nextflow.io/)
[Docker](https://www.docker.com/)

As all of the available parameters can already be found within the [`nextflow.config`](nextflow.config) file, an example run of the pipeline on the command line may look like this:

```bash
nextflow run lifebit-ai/pgx-nf
```

You can change the input parameters by either changing them in the [`nextflow.config`](nextflow.config) file or specifying them on the command line eg:
```bash
nextflow run main.nf --celllines testdata/annotation/cell_annotation_all.csv
```

## Deploit

Deploit is a bioinformatics platform, developed by Lifebit, where you can run your analysis over the Cloud/AWS.

It is free for indivudal users to [sign-up](https://deploit.lifebit.ai/register)

To run the pipeline once logged in navigate to the pipelines page:

![deploit](https://raw.githubusercontent.com/lifebit-ai/ecw-converter/master/images/deploit.png)

### Import the pipeline
To import the pipeline into Deploit you simply need to copy & paste the GitHub repository URL: `https://github.com/lifebit-ai/pgx-nf`. To do so follow the steps below:
![import_pgx_deploit](https://raw.githubusercontent.com/lifebit-ai/images/master/pgx-nf/import_pgx.gif)


### Running on Deploit

![run_pgx_deploit](https://github.com/lifebit-ai/images/blob/master/pgx-nf/run_pgx.gif)

The pipeline can then be run in three simple steps:
1. **Find & select the pipeline** After importing the pipeline into Deploit you are immediately taken to the input data & parameters page and so might not need to do this. Imported pipelines are found under the `MY PIPELINES & TOOLS` section. The name with be the same as you used to import the pipeline, eg `pgx-nf`. 
2. **Enter the data & parameters** here you can specify the input parameters & data. You may not need specify any parameters as they are already provided in the [`nextflow.config`](nextflow.config) file. However, you can overwrite the parameters in the config file by specifying parameters on Deploit
3. **Select a project & instance** before running the analysis you must select the instance. This determines the available resources & cost. The project allows multiple jobs to be grouped together/organised

See [example run](https://deploit.lifebit.ai/public/jobs/5d5687c7266ca500c4085c9b) of the pipeline on Deploit:
[![pgx_job_page](https://raw.githubusercontent.com/lifebit-ai/images/master/pgx-nf/pgx_job_page.png)](https://deploit.lifebit.ai/public/jobs/5cc2e65702877100b2bb96e6)