rule sno_host_classification:
    input:
        tableS1 = join(config['path']['tmp'],
                       config['file']['tableS1']),
        snodb = join(config['path']['ref'],
                     config['file']['snoDB']),
        bio_function = join(config['path']['ref'],
                            config['file']['bio_function_completed']),
        gab_sno_host = join(config['path']['tmp'],
                            config['file']['gab_sno_host']),
    output:
        'tok'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/sno_host_classification.py"
