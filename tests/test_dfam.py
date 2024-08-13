from seqbank.dfam import dfam_url

def test_dfam_url_default():
    """Test the default parameters."""
    expected_url = "https://www.dfam.org/releases/current/families/Dfam_curatedonly.h5.gz"
    assert dfam_url() == expected_url


def test_dfam_url_curated_true():
    """Test when curated is True (explicitly set) and release is 'current'."""
    expected_url = "https://www.dfam.org/releases/current/families/Dfam_curatedonly.h5.gz"
    assert dfam_url(curated=True, release="current") == expected_url


def test_dfam_url_curated_false():
    """Test when curated is False."""
    expected_url = "https://www.dfam.org/releases/current/families/Dfam.h5.gz"
    assert dfam_url(curated=False, release="current") == expected_url

def test_dfam_url_custom_release():
    """Test with a custom release value."""
    expected_url = "https://www.dfam.org/releases/2023-07/families/Dfam_curatedonly.h5.gz"
    assert dfam_url(curated=True, release="2023-07") == expected_url

def test_dfam_url_curated_false_custom_release():
    """Test when curated is False with a custom release value."""
    expected_url = "https://www.dfam.org/releases/2023-07/families/Dfam.h5.gz"
    assert dfam_url(curated=False, release="2023-07") == expected_url
