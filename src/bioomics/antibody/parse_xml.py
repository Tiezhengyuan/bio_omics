import os
import xml.etree.ElementTree as ET
import re


class ParseXml:

    @staticmethod
    def href_text(xml_string):
        res = []
        xml_string = xml_string.replace('&', '')
        root = ET.fromstring(xml_string)
        for a in root.findall('a'):
            rec = {}
            try:
                rec['text'] = root.text
                rec['value'] = a.text
                rec['url'] = a.get('href')
            except Exception as e:
                pass
            res.append(rec)
        return res

    @staticmethod
    def ul_li(xml_string):
        res = []
        xml_string =xml_string.replace('""', '"')
        xml_string =xml_string.replace('<a>', '</a>')
        xml_string = re.sub(r'&|®|<i>|</i>|<br>|</br>', '', xml_string)
        try:
            root = ET.fromstring(xml_string)
            for ul in root.findall('ul'):
                for li in ul.findall('li'):
                    rec = {
                        'text': li.text,
                    }
                    try: 
                        rec['value'] = li.find('a').text
                        rec['url'] = li.find('a').get('href')
                    except Exception as e:
                        pass
                    res.append(rec)
        except Exception as e:
            pass
        return res

    @staticmethod
    def tag_substr(text, head_str, end_tag):
        if head_str in text:
            start = text.index(head_str)
            substr = text[start+len(head_str):]
            if end_tag in substr:
                end = substr.index(end_tag)
                substr = substr[:end+len(end_tag)]
                return substr.replace('\n', '')