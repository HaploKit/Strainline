
import sys
import os

def compute_ANI_ij(fa,i,j,outdir):
    '''
    compute the ANI of the i <-> j th sequence in the fasta (one sequence per line)
    '''
    try:
        i_fa="{}/seq.{}.fa".format(outdir,i)
        j_fa="{}/seq.{}.fa".format(outdir,j)
        os.system("head -{} {}|tail -2 >{}".format(i*2,fa,i_fa))

        os.system("head -{} {}|tail -2 >{}".format(j*2,fa,j_fa))
        os.system("fastANI -q {} -r {} -o {}/ani.{}.{}.out -t 1 >/dev/null 2>&1".format(i_fa,j_fa,outdir,i,j))
        os.system("rm -f {} {}".format(i_fa,j_fa))
    except:
        raise Exception("Failed to extract sequences for ANI computation.")


if __name__ == '__main__':
    fa,i,j,outdir=sys.argv[1:]
    i =int(i)
    j =int(j)
    os.system("mkdir -p {}".format(outdir))
    compute_ANI_ij(fa,i,j,outdir)
