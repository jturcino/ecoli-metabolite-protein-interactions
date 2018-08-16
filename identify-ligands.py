#!/usr/bin/env python3

# setup
from argparse import ArgumentParser
from zipfile import ZipFile
from os.path import isfile
from re import split, search
from numpy import empty, nan
from xml.etree import ElementTree
from requests import get
import pandas as pd

# global variables
ecocyc = 0
hmdb = 1
invalid_mids = ['NA', 'multiple charge', 'neg', 'nan']
skipsheets = ['Samplelist of essential', 'TFs']

# HMDB lookup done via XML
hmdbfile = 'hmdb_metabolites.zip'
print('Loading HMDB from', hmdbfile)
assert isfile(hmdbfile), 'HMDB file ({}) does not exist'.format(hmdbfile)
zf = ZipFile(hmdbfile, 'r')
f = zf.open('{}.xml'.format(hmdbfile[:-4]))
hmdbroot = ElementTree.parse(f).getroot()
zf.close()

def hmdb_lookup(mid, xmlns='{http://www.hmdb.ca}'):
    '''Uses HMDB XML to extract InChI, InChIKey, and SMILES for a given metabolite ID'''
    smile = inchi = inchikey = None
    e = hmdbroot.find("./{}metabolite/[{}name='{}']".format(xmlns, xmlns, mid))
    if e is not None:
        smile = e[16].text
        inchi = e[17].text[6:] # stripping 'InChI=' from front
        inchikey = e[18].text
    return dict(InChI=inchi, InChIKey=inchikey, SMILES=smile, MetaboliteID=mid)

# EcoCyc  done via cURL
def ecocyc_lookup(mid, org='ECOLI', url='https://websvc.biocyc.org/getxml?'):
    '''Uses requests.get and regex to extract InChI, InChIKey, and SMILES for a given metabolite ID'''
    # get info
    url += '{}:{}'.format(org, mid)
    resp = get(url)
    resp.raise_for_status()
    resp = resp.text
    # pull values from XML values (None if not found)
    inchi = search(r'(?<=InChI=)[^<]*(?=</inchi>)', resp)
    inchikey = search(r'(?<=InChIKey=).*(?=</inchi-key>)', resp)
    smile = search(r'(?<=<string title=\'smiles\'>).*(?=</string>)', resp)
    # if all values found (not None), parse
    if None not in [inchi, inchikey, smile]:
        inchi = inchi.group()
        inchikey = inchikey.group()
        smile = smile.group()
    return dict(InChI=inchi, InChIKey=inchikey, SMILES=smile, MetaboliteID=mid)

# general functions
def init_file_dfs(subheaders):
    '''Inits pandas.DataFrame for the nodata, simple, match, multimatch, and ambiguous output Excel file sheets. Uses subheaders pandas.Series to populate column names and first row'''
    # nodata.xlsx for rows without valid EcoCyc or HMDB entries or no matches between databases
    # drop EcoCyc and HMDB columns
    nodata = subheaders.drop(['EcoCyc', 'HMDB']).to_frame().transpose()
    
    # match.xlsx for rows with single overlapping valid value between databases
    # copy nodata
    match = nodata.copy()
    # init InChI, InChIKey, and SMILES columns with NaN
    match['InChI'] = match['InChIKey'] = nan
    # init MetaboliteID columns with respective database names
    match['eMetaboliteID'] = match['eSMILES'] = 'EcoCyc'
    match['hMetaboliteID'] = match['hSMILES'] = 'HMDB'

    # multimatch.xlsx for rows with multiple matches
    # copy match
    # multiple values will be separated with "|"
    multimatch = match.copy()
    
    # simple.xlsx for rows with only one database and one valid value
    # copy nodata 
    simple = nodata.copy()
    # init InChI, InChIKey, SMILES, MetaboliteID, and Database columns with NaN
    simple['InChI'] = simple['InChIKey'] = simple['SMILES'] = simple['MetaboliteID'] = simple['Database'] = nan
    
    # ambiguous.xlsx for rows with no matches and multiple values
    # copy nodata 
    ambiguous = nodata.copy()
    # init InChI, InChIKey, SMILES, and MetaboliteID with database names
    ambiguous['eMetaboliteID'] = ambiguous['eInChI'] = ambiguous['eInChIKey'] = ambiguous['eSMILES'] = 'EcoCyc'
    ambiguous['hMetaboliteID'] = ambiguous['hInChI'] = ambiguous['hInChIKey'] = ambiguous['hSMILES'] = 'HMDB'

    # return dataframes
    return nodata, simple, match, multimatch, ambiguous

def getlists(excelrow):
    '''Takes row of excel data file and returns lists of EcoCyc and HMDB metabolite IDs'''
    assert 'HMDB' in excelrow.keys() and 'EcoCyc' in excelrow.keys(), 'DB columns not in provided row'
    ids = {}
    for i in ['EcoCyc', 'HMDB']:
        # extract metabolite strings and split
        ids[i] = str(excelrow[i])
        ids[i] = split(r"; |/(?![\w:/\-]+\))", ids[i])
        # throw out invalid values
        ids[i] = [ j for j in ids[i][j] if j not in invalid_mids and not j.isdigit() ]
    # return lists of IDs
    return ids['EcoCyc'], ids['HMDB']

def retrieve(mids, db, columns=['InChI', 'MetaboliteID', 'InChIKey', 'SMILES']):
    '''Retrieves InChI, InChIKey, and SMILES from a list of metabolite IDs via XML database lookup. Returns pandas.DataFrame'''
    assert db == hmdb or db == ecocyc, 'Invalid database: {}'.format(db)
    # init data frame with rows for each mid in mids and four columns
    df = pd.DataFrame(empty([len(mids),len(columns)], dtype=str), columns=columns, dtype=str)
    # fill each row
    for i in range(len(mids)):
        df.iloc[i] = (ecocyc_lookup(mids[i]) if db == ecocyc else hmdb_lookup(mids[i]))
    # return with empty values removed
    return df.dropna()

def getmatches(einfo, hinfo, columns=['InChI', 'InChIKey', 'hSMILES', 'eSMILES', 'eMetaboliteID', 'hMetaboliteID']):
    '''Parses matching metabolites using InChIKeys. Takes pandas.DataFrame for each database and returns pandas.DataFrame if match(es) found. If no matches found, returns None'''
    # find intersecting InChIKeys
    matchkeys = list(set(einfo['InChIKey']) & set(hinfo['InChIKey']))
    # init return value and fill if matches exist
    df = (None if len(matchkeys)==0 else pd.DataFrame(empty([len(matchkeys), len(columns)], dtype=str), columns=columns, dtype=str))
    for i in range(len(matchkeys)):
        key = matchkeys[i]
        # assign values from hinfo (InChI and InChIKey will be same)
        d = hinfo.loc[hinfo['InChIKey']==key].iloc[0].to_dict()
        d['hMetaboliteID'] = d.pop('MetaboliteID')
        d['hSMILES'] = d.pop('SMILES')
        # assign values from einfo
        d['eMetaboliteID'] = einfo.loc[einfo['InChIKey']==key]['MetaboliteID'][0]
        d['eSMILES'] = einfo.loc[einfo['InChIKey']==key]['SMILES'][0]
        # add values to columns
        df.iloc[i] = d
    return df

# row building functions
def joinvalues(x, sep='|'):
    '''Joins values in pandas.Series with sep. Returns string.'''
    return sep.join(list(x))

def process_ogrow(ogrow):
    '''Removes EcoCyc and HMDB columns from original excel sheet row and returns remaining values as dict.'''
    ogrow.drop(['EcoCyc', 'HMDB'], inplace=True)
    return ogrow.to_dict()

def build_simplerow(ogrow, minfo, db):
    '''Takes row from original excel sheet and returns row for simple.xlsx'''
    # remove EcoCyc and HMDB entries from row
    # then convert to dict
    d = process_ogrow(ogrow)
    # add columns and new values
    d.update(minfo.iloc[0].to_dict())
    d['Database'] = db
    return d

def getdict_ambiguousinfo(ambig_df, db, sep='|'):
    '''Helper function for formatting ambiguous data for file rows. Joins ambiguous values with sep parameter. Returns dict.'''
    assert db==ecocyc or db==hmdb, 'Invalid database: {}'.format(db)
    # add e or h to column names to indicate db
    prefix = ('e' if db==ecocyc else 'h')
    ambig_df.columns = [prefix+i for i in ambig_df.columns]
    # build dict and return
    return { i : joinvalues(ambig_df[i]) for i in ambig_df }

def build_ambiguousrow(ogrow, minfo1, db1, minfo2=None, db2=None):
    '''Takes row from original excel sheet and returns row for ambiguous.xlsx. Second database set optional.'''
    # remove EcoCyc and HMDB entries from row
    # then convert to dict
    d = process_ogrow(ogrow)
    # add columns and new values for first set
    # repeat with second set if provided
    d.update(getdict_ambiguousinfo(minfo1, db1))
    if minfo2 is not None and db2 is not None:
        d.update(getdict_ambiguousinfo(minfo2, db2))
    return d

def build_matchrow(ogrow, match):
    '''Takes row from original excel sheet and returns row for match.xlsx. or matches.xlsx'''
    # remove EcoCyc and HMDB entries from row
    # then convert to dict
    d = process_ogrow(ogrow)
    # add columns and new values to dict
    d.update({ i : joinvalues(match[i]) for i in match })
    return d

if __name__ == '__main__':

    # set up args
    parser = ArgumentParser()
    parser.add_argument('-e', '--excelfile', dest='excelfile', default='raw-data/essential_2017Nov(July2018).xlsx', help='Raw data excel file')
    parser.add_argument('-o', '--outdir', dest='outdir', default='processed-data', help='Output directory')
    args = parser.parse_args()

    # init dictionary of dataframes to be written as output files
    # as each sheet of the original excel file is processed, add (sheetname, dataframe) tuples to list
    fdict = { 'nodata':[], 'simple':[], 'match':[], 'multimatch':[], 'ambiguous':[] }

    # load excel file
    print('Loading raw file ({})'.format(args.excelfile))
    data = pd.ExcelFile(args.excelfile)
    tabs = data.sheet_names

    # process each sheet in data excel file
    for sheetname in tabs:
        # skip sheets not containing metabolite info; otherwise, load sheet data
        if sheetname in skipsheets:
            print('Skipping', sheetname)
            continue
        print('Processing', sheetname)
        sheet = data.parse(sheetname=sheetname, skiprows=0)

        # drop 'Unnamed' columns and init dataframes to be saved in fdict before writing
        cols = [i for i in sheet.columns if i[:7]!='Unnamed']
        sheet = sheet[cols]
        assert 'EcoCyc' in sheet.columns and 'HMDB' in sheet.columns, 'Sheet {} lacks EcoCyc and HMDB columns:\n{}'.format(sheetname, sheet.columns)
        subheaders = sheet.iloc[0]
        nodata_df, simple_df, match_df, multimatch_df, ambiguous_df = init_file_dfs(subheaders)

        # process rows
        for i in range(1,len(sheet)):
            row = sheet.iloc[i]
            # get lists of EcoCyc and HMDB metabolite IDs
            emids, hmids = getlists(row)
            
            # skip rows lacking valid entries
            if len(emids)==0 and len(hmids)==0: 
                print('\tSkipping row', i)
                rowdict = process_ogrow(row)
                nodata_df = nodata_df.append(rowdict, ignore_index=True)
                continue

            # if only one database is empty, get info
            elif len(emids)==0 or len(hmids)==0:
                # set id list and database name; then get info
                ids, database = (hmids, hmdb if len(emids)==0 else emids, ecocyc)
                info = retrieve(ids, database)
                # if only one valid value, store in simple.xlsx; else store in ambiguous.xlsx
                if info.shape[0] == 1:
                    print('\tSimple row', i)
                    rowdict = build_simplerow(row, info, database)
                    simple_df = simple_df.append(rowdict, ignore_index=True)
                else:
                    print('\tAmbiguous row', i, '(single database)')
                    rowdict = build_ambiguousrow(row, info, database)
                    ambiguous_df = ambiguous_df.append(rowdict, ignore_index=True)

            # if both databases have entries, get info for both
            else:
                # get info on metabolite IDs
                ecocycinfo = retrieve(emids, ecocyc)
                hmdbinfo = retrieve(hmids, hmdb)
                # find matches
                matches = getmatches(ecocycinfo, hmdbinfo)
                # if match(es), add to (multi)match_df
                if matches is not None:
                    rowdict = build_matchrow(row, matches)
                    if len(matches)==1:
                        print('\tSingle match row', i)
                        match_df = match_df.append(rowdict, ignore_index=True)
                    else:
                        print('\tMultiple match row', i)
                        multimatch_df = multimatch_df.append(rowdict, ignore_index=True)
                # whether or not any matches, add to ambiguous_df
                rowdict = build_ambiguousrow(row, ecocycinfo, ecocyc, hmdbinfo, hmdb)
                ambiguous_df = ambiguous_df.append(rowdict, ignore_index=True)
        # END PROCESS ROWS

        # now done building dfs, so append to fdict value lists
        # format as tuples (sheetname, df)
        print('Appending', sheetname, 'dataframes')
        fdict['nodata'].append((sheetname, nodata_df))
        fdict['simple'].append((sheetname, simple_df))
        fdict['match'].append((sheetname, match_df))
        fdict['multimatch'].append((sheetname, multimatch_df))
        fdict['ambiguous'].append((sheetname, ambiguous_df))
    # END PROCESS SHEETS

    # write files
    print('Writing files ({})'.format(args.outdir))
    for name, tuplelist in fdict.items():
        fname = '{}/{}.xlsx'.format(args.outdir, name)
        print('\tWriting', fname)
        writer = pd.ExcelWriter(fname, engine='xlsxwriter')
        for sheetname, df in tuplelist:
            df.to_excel(writer, sheet_name=sheetname, index=False)
        writer.save()
    # END WRITE FILES