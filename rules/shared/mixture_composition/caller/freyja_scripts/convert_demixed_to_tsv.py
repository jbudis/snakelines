import sys
import pandas as pd

def convert_to_tsv(file_path, sample_name, out_path):
    occurrences = {}
    sample_res = {}
    line = ''
    try:
        f = open(file_path)
        line = f.read()
        f.close()
    except:
        sys.exit(1)
    
    lineages = line.split('\n')[2].split()[1:]
    abundances = line.split('\n')[3].split()[1:]
    assert len(lineages) == len(abundances)

    for i in range(len(lineages)):
        sample_res[lineages[i]] = float(abundances[i])
    
    occurrences[sample_name] = sample_res

    occ_df = pd.DataFrame.from_dict(occurrences, orient='index')
    occ_df.to_csv(out_path)

if __name__ == "__main__":
    convert_to_tsv(sys.argv[1], sys.argv[2], sys.argv[3])
