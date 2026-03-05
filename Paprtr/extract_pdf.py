import sys
import pypdf

reader = pypdf.PdfReader(sys.argv[1])
text = chr(10).join([page.extract_text() for page in reader.pages if page.extract_text()])

with open(sys.argv[2], 'w', encoding='utf-8') as f:
    f.write(text)
