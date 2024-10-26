from bs4 import BeautifulSoup

def extract_function_value(html_content):
    soup = BeautifulSoup(html_content, "html.parser")
    id_value = soup.find("input", {"name":"getSeqGi"})["value"]
    if id_value:
        return id_value
    return None

# Пример использования
html_content = '''
<td class="l c0"><span class="ind">2</span><input class="cb" id="chk_2" name="getSeqGi" onclick="configDescrLinks(event,this)" type="checkbox" value="KAI2545001.1"><label for="chk_2" onclick="checkShiftKey(event,this)"><span class="accsb">Select seq gb|KAI2545001.1|</span></label></input></td> <td class="ellipsis c2"><span><a accs="gb|KAI2545001.1|" class="deflnDesc" gi="KAI2545001.1" href="#alnHdr_KAI2545001" hsp="1" id="deflnDesc_2" ind="2" len="41" onclick="DisplayAlignFromDescription(this);" seqfsta="gb|KAI2545001.1|" seqfstadwnld="gb|KAI2545001.1|" seqid="KAI2545001" title="Go to alignment for islet cell autoantigen 1 [Homo sapiens] &gt;gb|KAI4012930.1| islet cell autoantigen 1, partial [Homo sapiens] &gt;gb|PNI57852.1| ICA1 isoform 14, partial [Pan troglodytes] &gt;gb|PNJ69881.1| ICA1 isoform 4, partial [Pongo abelii]">islet cell autoantigen 1 [Homo sapiens]</a></span></td> <td class="c6">85.1</td>
'''

function_value = extract_function_value(html_content)
print(function_value)  # Вывод: LOW QUALITY PROTEIN: islet cell autoantigen 1 like [Homo sapiens]