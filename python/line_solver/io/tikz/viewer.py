"""
Opens PDF files generated from TikZ diagrams in an external viewer.

Port of jar/src/main/java/jline/io/tikz/TikZViewer.java. Standalone viewers
are preferred over the desktop default, which may open a browser.
"""

import os
import shutil
import subprocess
import sys

_STANDALONE_VIEWERS = ('evince', 'okular', 'mupdf', 'zathura', 'qpdfview', 'acroread')


def _desktop_open_command():
    """The platform's open-with-default-application command, or None."""
    if sys.platform == 'darwin':
        return ['open']
    if sys.platform.startswith('win'):
        return None
    return ['xdg-open'] if shutil.which('xdg-open') else None


def display_pdf(pdf_file, title=None):
    """Open pdf_file in a standalone viewer, else the desktop default; print the path when neither works."""
    pdf_file = os.path.abspath(pdf_file)
    if not os.path.exists(pdf_file):
        print('PDF file does not exist: ' + pdf_file, file=sys.stderr)
        return False
    for viewer in _STANDALONE_VIEWERS:
        if shutil.which(viewer):
            try:
                subprocess.Popen([viewer, pdf_file], stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                                 stderr=subprocess.DEVNULL, start_new_session=True)
                if title is not None:
                    print('Opened diagram for: ' + title)
                return True
            except OSError:
                continue
    try:
        if sys.platform.startswith('win'):
            os.startfile(pdf_file)  # noqa: the Windows desktop default
            opened = True
        else:
            cmd = _desktop_open_command()
            opened = False
            if cmd is not None:
                subprocess.Popen(cmd + [pdf_file], stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                                 stderr=subprocess.DEVNULL, start_new_session=True)
                opened = True
        if opened:
            if title is not None:
                print('Opened diagram for: ' + title)
            return True
    except OSError as exc:
        print('Failed to open PDF with default viewer: ' + str(exc), file=sys.stderr)
    print('Could not open PDF viewer. File saved at: ' + pdf_file)
    return False


def is_pdf_viewer_available():
    """True when a standalone viewer or a desktop open command exists."""
    if any(shutil.which(v) for v in _STANDALONE_VIEWERS):
        return True
    return sys.platform.startswith('win') or _desktop_open_command() is not None


displayPDF = display_pdf
isPdfViewerAvailable = is_pdf_viewer_available
