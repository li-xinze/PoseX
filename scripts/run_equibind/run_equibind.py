import os
import argparse
import yaml


def main(args: argparse.Namespace):
    os.chdir(os.path.join(args.equibind_exec_dir))
    with open(args.yml_path, 'r', encoding='utf-8') as f:
        file_content = f.read()
    data = yaml.load(file_content, yaml.FullLoader)
    data["inference_path"] = args.input_dir
    data["output_directory"] = args.output_dir

    modified_yml = os.path.join(os.path.dirname(args.yml_path), "modified_inference.yml")
    with open(modified_yml, 'w', encoding='utf-8') as f:
        yaml.dump(data, stream=f, allow_unicode=True, encoding='utf-8')
    cmd = f"python inference.py --config={modified_yml}"
    os.system(cmd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", type=str, required=True, help="Path to the input files")
    parser.add_argument("--output_dir", type=str, required=True, help="Path to save the output")
    parser.add_argument("--equibind_exec_dir", type=str, required=True, help="Path to the Equibind codebase")
    parser.add_argument("--yml_path", type=str, required=True, help="Path to the yml config file")
    args = parser.parse_args()

    main(args)
