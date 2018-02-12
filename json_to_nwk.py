#########################################################

# James Stimson released this initial version 12/02/2018
# Turns JSON tree into newick tree format

#########################################################

import sys

if len(sys.argv) == 1:
    print 'Please supply .json file name as argument'
    sys.exit()

jsonFileName = sys.argv[1]
jsonFile = open(jsonFileName,'r')
nwkFileName = jsonFileName[:-4] + 'nwk'
nwkFile = open(nwkFileName,'w')

depth = 0
brack_depth = 0
node_depth = 0
id_count = 0
t_count = 0
clade_count = 0
date_count = 0
raw_date_count = 0
brack_count = 0
end_brack_count = 0
nwkstr = ''
tvalue_used = False

remove_end_brace = False
t_depth = 0

node_times = []
node_levels = []
node_lengths = []
node_type = []

for line in jsonFile:
    if '{' in line:
        depth += 1
        if (line.rstrip()).lstrip() == '{':
                        
            if (len(nwkstr)>1) and (nwkstr[-1] != '('):
                nwkstr += ','
                 
            nwkstr += '('
            
            brack_count += 1
            brack_depth += 1
            
            
    if '}' in line:
        depth -= 1
        while brack_depth > depth:

            if remove_end_brace:
                remove_end_brace = False
            else:
                nwkstr += ')'
                brack_count -= 1
                end_brack_count -= 1

            if not tvalue_used:
                nwkstr += ':' + tvalue
                tvalue_used = True
                node_times.append(float(tvalue))
                node_levels.append(t_depth + 1)
                node_type.append('internal')
            
            end_brack_count += 1
            brack_depth -= 1
    if '"clade"' in line:
        rhs = (line.split(':')[1]).lstrip()
        clade_count += 1
    if '"strain"' in line:
        rhs = (line.split(':')[1]).lstrip()
        node_name = rhs.split('"')[1]
    if 'isolate_id' in line:
        id_count += 1
        rhs = line.split(':')[1]
        nwkstr = nwkstr[:-1] # could check this is indeed a bracket
        nwkstr += rhs.split('"')[1]
    if '"raw_date"' in line:
        raw_date_count += 1
        # rely on the fact that this comes shortly after isolate_id
        rhs = line.split(':')[1]    
        if tvalue_used:
            print 'value used error'
        nwkstr += '_' + rhs.split('"')[1] + ':' + tvalue
        node_times.append(float(tvalue))
        node_levels.append(depth)
        node_type.append('tip')
        remove_end_brace = True
        tvalue_used = True
        node_depth = depth
    if 'tvalue' in line:
        rhs = (line.split(':')[1]).lstrip()
        tvalue = rhs.split(',')[0]
        tvalue_used = False
        t_depth = depth
        t_count += 1
    if 'num_date' in line:
        date_count += 1
    

# Infer lengths from times
for n in range(len(node_times)):
    attempt = n+1
    node_lengths.append(0.0)
    while attempt < len(node_times):
        if (node_type[attempt] == 'internal') and (node_levels[n] > node_levels[attempt]):
            node_lengths[n] = node_times[n] - node_times[attempt]
            break  
        attempt += 1
    
explode = nwkstr.split(':')
newnwkstr = explode[0]

for n in range(len(node_lengths)):
    newnwkstr += ':' + str(node_lengths[n]) + (explode[n+1])[len(str(node_times[n])):]
    
nwkFile.write(newnwkstr)
nwkFile.write(':0.0;')

print 'found ' + str(id_count) + ' isolates'
print 'found ' + str(raw_date_count) + ' raw dates'
print 'found ' + str(t_count) + ' times'
print 'found ' + str(date_count/2) + ' dates'
print 'found ' + str(clade_count) + ' clades'
#print 'created ' + str(brack_count) + ' of these ('
#print 'created ' + str(end_brack_count) + ' of these )'
#print 'this should be zero: ' + str(depth)

jsonFile.close()
nwkFile.close()
