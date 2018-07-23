%getGeneTable will get IMGT gene tables listing the gene name and host
%strain that contained the gene. The website can be found at 
%http://www.imgt.org/IMGTrepertoire/LocusGenes/#F
function getGeneTable()
disp('Code is under construction.')
return
URL = 'http://www.imgt.org/IMGTrepertoire/index.php?section=LocusGenes&repertoire=genetable&species=Mus_musculus&group=IGHV';

PageData = webread(URL);

TableStart = regexpi(PageData,'<table>','start');
TableEnd = regexpi(PageData,'</table>','end');
TableData = PageData(TableStart:TableEnd);

%For this table, reduce to the follow format

%Count rows


HeaderS = regexp(TableData,'<th[^\>]*>','end') + 1;
HeaderE = regexp(TableData,'</th>','end') + 1;
RowS = regexp(TableData,'<tr[^\>]*>','end') + 1;
RowE = regexp(TableData,'</tr>','start') - 1;
CellS = regexp(TableData,'<td[^\>]*>','end') + 1;
CellE = regexp(TableData,'</td>','start') - 1;




CellCount = regexp(TableData,'</td>','start');
Data = cell(length(RowCount) + length(CellCount),2);
