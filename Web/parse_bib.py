def count_substr(s, ss):
    c = 0
    for e in s:
        if e == ss:
            c += 1
    return c

class BibTeXerClass:
    
    def __init__(self, fName = '../CoolProp/CoolPropBibTeXLibrary.bib'):
        self.parse(fName)
        
    def _entries(self, lines):
        entries = []
        i = 0
        while i < len(lines):
            line = lines[i]
            entry = []
            if line.startswith('@') and not line.startswith('@comment'):
                iEnd = i
                openbrackets = count_substr(lines[iEnd],'{')-count_substr(lines[iEnd],'}')
                while openbrackets > 0:
                    openbrackets += count_substr(lines[iEnd+1],'{')-count_substr(lines[iEnd+1],'}')
                    iEnd += 1
                
                entry = ''.join(lines[i:iEnd+1]).replace('\n','').replace('\t',' ')
                head,kwargs = entry.split(',',1)
                d = {}
                
                d['family'], d['key'] = head.replace('@','').split('{')
                
                kwargs  = kwargs.rstrip('}')
                kwargs = kwargs.split('},')
                
                for kwarg in kwargs:
                    k,v = kwarg.split('= {')
                    
                    d[k.strip()] = v.replace('{','').replace('}','')
                    
                entries.append(d)
                
                i = iEnd+1
            else:
                i += 1
        return entries
                
    def entry2rst(self, key):
        
        entry = self.findentry(key)
        
        if entry is None:
            return ''
        
        if entry['family'] == 'ARTICLE':
            if 'journal' not in entry: entry['journal'] = ''
            if 'volume' not in entry: entry['volume'] = ''
            if 'pages' not in entry: entry['pages'] = ''
            
            return entry['author'] + ', ' + entry['year'] + ', ' + entry['title'] + ', *' + entry['journal'] + '*, ' + entry['volume'] + ':' + entry['pages']
        elif entry['family'] == 'CONFERENCE':
            if 'journal' not in entry: entry['journal'] = ''
            return entry['author'] + ', ' + entry['year'] + ', ' + entry['title'] + ', *' + entry['journal'] + '*'
        elif entry['family'] == 'MASTERSTHESIS':
            return entry['author'] + ', ' + entry['year'] + ', ' + entry['title'] + ', *' + entry['school'] + '*'
        elif entry['family'] == 'UNPUBLISHED':
            return entry['author'] + ', ' + entry['year'] + ', ' + entry['title'] + ', note: ' + entry['note']
        elif entry['family'] == 'BOOK':
            return entry['author'] + ', ' + entry['year'] + ', *' + entry['title'] + '*, ' + entry['publisher']
        elif entry['family'] == 'TECHREPORT':
            return entry['author'] + ', ' + entry['year'] + ', *' + entry['title'] + '*, ' + entry['institution']
        else:
            print entry
            raise ValueError(entry['family'])
        
    def findentry(self, key):
        for entry in self.entries:
            if entry['key'] == key:
                return entry
            
    def parse(self, fName):
        lines = open(fName,'r').readlines()
        self.entries = self._entries(lines)

if __name__=='__main__':
    B = BibTeXerClass()
    print B.entry2rst('Assael-JPCRD-2013A')