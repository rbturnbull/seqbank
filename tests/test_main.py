import tempfile
from typer.testing import CliRunner

from seqbank.main import app
from pathlib import Path

runner = CliRunner()


test_data = Path(__file__).parent / "testdata"


def test_export():
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmpdirname = Path(tmpdirname)

        exported_path = tmpdirname/"test_main_export.fasta"
        result = runner.invoke(app, ["export", str(test_data/"seqbank.sb"), str(exported_path)])
        assert result.exit_code == 0
        assert exported_path.exists()
        text = exported_path.read_text()
        assert text.count(">") == 6
        assert ">NC_024664.1\nTAAAAAGAAAAAAATTTT" in text
        assert ">NZ_JAJNFP010000161.1\nTAATATTTGCTTTCATTTCTAAATAG" in text
        assert ">NC_010663.1\nAGTTTTAAACCTCTGATCGAAC" in text
        assert ">NC_036112.1\nGAAATACCCAATATC" in text
        assert ">NC_036113.1\nGAAATACCCAATATCTTGTTCC" in text
        assert ">NC_044840.1\nGGGCGAACGACGGGAATTGAACCCGC" in text
        assert len(text) > 5000
