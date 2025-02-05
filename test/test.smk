def return_fils():
    return ['result/x_a1_b1.out',
            'result/x_a2_b1.out',
            'result/y_a2_b2.out']

rule all:
    input:
        return_fils()

rule a1:
    input:
        'input/{sample}.in'
    output:
        'tmp/{sample}_a1.out'
    shell:
        'echo done > {output}'

rule a2:
    input:
        'input/{sample}.in'
    output:
        'tmp/{sample}_a2.out'
    shell:
        'echo done > {output}'

rule b1:
    input:
        'tmp/{sample}_{a_choice}.out'
    output:
        'result/{sample}_{a_choice}_b1.out'
    shell:
        """
        cat {input} <(echo "done2") > {output}
        """

rule b2:
    input:
        'tmp/{sample}_{a_choice}.out'
    output:
        'result/{sample}_{a_choice}_b2.out'
    shell:
        """
        cat {input} <(echo "done2") > {output}
        """
