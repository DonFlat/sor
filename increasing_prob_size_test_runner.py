import argparse
import subprocess

# def generate_powers_of_two(n):
#     return [2**i for i in range(n + 1)]

def run_command(command):
    try:
        subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e.stderr}")

def main():
    parser = argparse.ArgumentParser(description="Process some named arguments.")

    parser.add_argument('--env', type=str, help='The environment to use, das or local')
    parser.add_argument('--app', type=str, help='raw | big | norm | split')
    parser.add_argument('-p', type=int, help='numbers of processor')

    args = parser.parse_args()

    env = args.env
    app = args.app
    node = args.p

    # [2^3, 2^25]
    sizes = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16834, 32768, 65536, 131072, 262144]
    if env == 'local':
        for size in sizes:
            local_cmd = f"mpirun -n {node} ./target/release/sor {app} {size} {node}"
            print(f'Running: {local_cmd}')
            run_command(local_cmd)
    elif env == 'das':
        for size in sizes:
            das_cmd = f"prun -np {node} -1 OMPI_OPTS=\"--mca btl tcp,self --mca btl_tcp_if_include ib0\" -script $PRUN_ETC/prun-openmpi `pwd`/./target/release/sor {app} {size} {node}"
            print(f'Running: {das_cmd}')
            run_command(das_cmd)
    else:
        print("Neither local | das")

if __name__ == "__main__":
    main()
