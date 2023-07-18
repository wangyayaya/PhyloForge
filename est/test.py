import os
import multiprocessing

def run_cmd(infile):
    cmd = f"for i in `grep '>' {infile} " \
          "|awk '{print $3}'|sort|uniq|sed 's/\[//g'|sed 's/\]//g'`;do " \
          f"grep -w $i {infile}" \
          "|head -n 1|awk '{print $1}'|sed 's/>//g' >>" \
          f"cds_uniq/`basename {infile} .cds`.uniq; done"
    return os.system(cmd)
    # print(cmd)


def run_cmd_mul():
    file_list = [f for f in os.listdir('cds_raw') if f != '']
    in_file_list = [os.path.join('cds_raw', infile) for infile in file_list]
    p = multiprocessing.Pool(19)
    p.map(run_cmd, in_file_list)
    p.close()
    p.join()


if __name__ == '__main__':
    run_cmd_mul()