"""Synthetic launcher tests only; no scientific script is executed."""
from pathlib import Path
import importlib.util
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.dont_write_bytecode=True
path=Path(__file__).resolve().parents[1]/'run.py'
spec=importlib.util.spec_from_file_location('release_launcher',path)
release=importlib.util.module_from_spec(spec); spec.loader.exec_module(release)

class LauncherTests(unittest.TestCase):
    def test_containment_rejects_escape(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(ValueError):release.contained(Path(d),'../escape')

    def test_default_is_check_only(self):
        with tempfile.TemporaryDirectory() as d:
            inp=Path(d)/'in'; out=Path(d)/'out'
            with patch.object(sys,'argv',['run.py','--input-root',str(inp),'--output-root',str(out)]), \
                 patch.object(release,'validate') as validate, patch.object(release.subprocess,'run') as execute:
                release.main();validate.assert_called_once();execute.assert_not_called()
            self.assertFalse(out.exists())

    def test_nested_output_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            with patch.object(sys,'argv',['run.py','--input-root',d,'--output-root',str(Path(d)/'out')]), \
                 patch.object(release,'validate') as validate:
                with self.assertRaises(SystemExit) as caught:release.main()
                self.assertEqual(caught.exception.code,2);validate.assert_not_called()

    def test_existing_stage_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'baseline').mkdir()
            with self.assertRaises(FileExistsError):release.command('baseline',root/'in',root,'Rscript')

    def test_missing_depth_inputs_do_not_create_output(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            with patch.object(release.shutil,'which',return_value='Rscript'):
                with self.assertRaises(FileNotFoundError):release.command('depth-statistics',root/'in',root,'Rscript')
            self.assertFalse((root/'depth_statistics').exists())

    def test_r_retry_uses_unchanged_definition(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            with patch.object(release.shutil,'which',return_value='Rscript'):
                cmd=release.command('daa-retry',root/'in',root,'Rscript')
            self.assertTrue(cmd[1].endswith('retry_ours_family_numeric.R'))
            self.assertTrue(cmd[5].endswith('daa_ours_family_retry_v1'))
            self.assertEqual(cmd[-2:],['--mode','daa'])

if __name__=='__main__':unittest.main()
