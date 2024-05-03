import sys

# Argumentos do script: [1] Arquivo de mapeamento, [2] Arquivo VCF original
map_file_path = sys.argv[1]
vcf_file_path = sys.argv[2]

# Ler o arquivo de mapeamento para construir o dicionário de renomeação
rename_map = {}
with open(map_file_path, 'r') as map_file:
    for line in map_file:
        old_name, new_name = line.strip().split('\t')
        rename_map[old_name] = new_name

# Gerar novo cabeçalho VCF
with open(vcf_file_path, 'r') as vcf:
    for line in vcf:
        if line.startswith('#CHROM'):
            columns = line.strip().split('\t')
            new_columns = [rename_map.get(col, col) for col in columns]
            print('\t'.join(new_columns))
        else:
            print(line.strip())
