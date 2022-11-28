import os
import types
import point_cloud_utils as pcu


api_reference_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                  "sections", "api_reference.md")
if os.path.exists(api_reference_path):
    os.remove(api_reference_path)

with open(api_reference_path, "w") as api_reference_f:
    api_reference_f.write("# API Reference\n\n")

    for attr_name in dir(pcu):
        if attr_name.startswith('_'):
            print(f"Skipping private function {attr_name}")
            continue
        if attr_name in ['warn', 'contextmanager', 'np', 'os']:
            print(f"Skipping builtin")
            continue

        if eval(f"type(pcu.{attr_name})") == types.ModuleType:
            print(f"Skipping module {attr_name}")
            continue
        # if eval(f"type(pcu.{attr_name})") not in (types.BuiltinFunctionType, types.:
        #     print(f"Skipping non builtin {attr_name}")
        #     continue

        docstring = eval(f"pcu.{attr_name}.__doc__")

        # print(docstring)
        if docstring is None:
            print(f"Skipping empty docstring {attr_name}")
            continue
        if len(docstring.strip()) == 0:
            print(f"Skipping empty docstring {attr_name}")
            continue

        min_indent = 1000000
        for line in docstring.split("\n"):
            if len(line) == 0:
                continue
            indent = len(line) - len(line.lstrip())
            min_indent = min(indent, min_indent)

        # line_lengths = [len(line) - len(line.lstrip()) for line in docstring.split("\n")]
        # min_indent = min(line_lengths)
        if min_indent > 0:
            new_docstring = ""
            for line in docstring.split("\n"):
                new_docstring += line[min_indent:] + "\n"
            docstring = new_docstring

        # Trim empty lines at the beginning
        new_docstring = ""
        found_nonempty = False
        for line in docstring.split("\n"):
            if len(line.strip()) > 0:
                found_nonempty = True

            if not found_nonempty and len(line.strip()) == 0:
                continue

            new_docstring += line + "\n"
        docstring = new_docstring

        # api_reference_f.write(f"#### {attr_name}\n")
        # api_reference_f.write("```text\n")
        # api_reference_f.write(docstring)
        # api_reference_f.write("```")
        # api_reference_f.write("\n\n")

        # api_reference_f.write(f"## {attr_name}\n")
        api_reference_f.write(f"::: point_cloud_utils.{attr_name}\n    handler: python\n    options:\n        show_root_heading: true\n        show_source: false\n\n")
        api_reference_f.write("\n\n----------------\n\n")



