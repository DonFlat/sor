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
    parser.add_argument('--app', type=str, help='raw | rma | norm')
    parser.add_argument('--size', type=int, help='problem size')

    args = parser.parse_args()

    env = args.env
    app = args.app
    size = args.size

    node_num = [2, 3, 4, 5, 6, 7, 8]

    if env == 'local':
        for node in node_num:
            local_cmd = f"mpirun -n {node} ./target/release/sor {app} {size} {node}"
            print(f'Running: {local_cmd}')
            run_command(local_cmd)
    elif env == 'das':
        for node in node_num:
            das_cmd = f"prun -np {node} -1 -script $PRUN_ETC/prun-openmpi `pwd`/./target/release/sor {app} {size} {node}"
            print(f'Running: {local_cmd}')
            run_command(das_cmd)
    else:
        print("Neither local | das")

if __name__ == "__main__":
    main()
