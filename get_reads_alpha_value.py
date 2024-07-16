import sys
import os
import subprocess


# def get_reads_in_block(wgbs_tool_path, block_file_path, pat_file, output_file):
#    wgbs = wgbs_tool_path + 'wgbstools cview --strict --min_len 4 -L ' + block_file_path
#    command = wgbs + ' ' + pat_file + ' -o ' + output_file
#    subprocess.run(command, shell=True, capture_output=True, text=True)


if __name__ == "__main__":

    # output directory
    wgbs_tool_path = sys.argv[1]
    block_file_path = sys.argv[2]
    pat_file = sys.argv[3]   ### pat.gz file
    output_file = sys.argv[4]
    
    #result = get_reads_in_block(wgbs_tool_path, block_file_path, pat_file, output_file)
    wgbs = wgbs_tool_path + '/wgbstools cview --strict --min_len 4 -L ' + block_file_path
    command = wgbs + ' ' + pat_file + ' -o ' + output_file
    subprocess.run(command, shell=True, capture_output=True, text=True)

    reads = []
    out_dir = os.path.dirname(output_file)
    name = os.path.basename(output_file).split('.bed')[0]
    outfile = out_dir + '/' + name + '_alpha.bed'
    with open(output_file, 'r') as f:
        for line in f:
            line = line.strip().split()
            seq = line[2]
            alpha = line[2].count('C')/len(line[2])
            line.append(str(alpha))
            new_line = '\t'.join(line)
            reads.append(new_line)
            reads.append('\n')
    with open(outfile, 'w') as out:
        out.writelines(reads)
    
    rm_bed = 'rm ' + output_file
    subprocess.run(rm_bed, shell=True, capture_output=True, text=True)


