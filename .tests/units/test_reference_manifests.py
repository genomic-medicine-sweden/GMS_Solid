import os
import re
import unittest

import yaml

CONFIG_DIR = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../config"))

# Data config and the reference manifests that download its files
DATA_CONFIGS = {
    "config.data.hg19.yaml": [
        "design_files.hg19.yaml",
        "references.hg19.yaml",
        "nextseq.hg19.pon.yaml",
        "novaseq.hg19.pon.yaml",
    ],
    "config.data.hg38.GMS-Solid.yaml": [
        "design_files.hg38.GMS-Solid.yaml",
        "references.hg38.GMS-Solid.yaml",
        "novaseq.hg38.GMS-Solid.pon.yaml",
    ],
}


def get_manifest_paths(entry):
    if isinstance(entry, dict):
        for key, value in entry.items():
            if key in ("path", "link"):
                yield value
            else:
                yield from get_manifest_paths(value)


class TestReferenceManifests(unittest.TestCase):
    def test_data_config_paths_in_manifests(self):
        for data_config, manifests in DATA_CONFIGS.items():
            with self.subTest(data_config=data_config):
                with open(os.path.join(CONFIG_DIR, data_config)) as f:
                    used = set(re.findall(r"\{\{PROJECT_[A-Z]+_DATA\}\}/([^\s\"',]+)", f.read()))
                provided = set()
                for manifest in manifests:
                    with open(os.path.join(CONFIG_DIR, "references", manifest)) as f:
                        provided.update(get_manifest_paths(yaml.safe_load(f)))
                # A path is also provided when it is inside a downloaded folder
                missing = sorted(
                    path for path in used
                    if path not in provided and not any(path.startswith(p.rstrip("/") + "/") for p in provided)
                )
                self.assertEqual(missing, [], f"{data_config} uses files not in {', '.join(manifests)}")


if __name__ == "__main__":
    unittest.main()
