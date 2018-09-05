#!/usr/bin/env python
#coding: utf-8

import pandas as pd
import ete3
import re
import os, sys

ncbi     = ete3.NCBITaxa()

header = 'assembly_accession bioproject biosample wgs_master refseq_category taxid species_taxid organism_name infraspecific_name isolate version_status assembly_level release_type genome_rep seq_rel_date asm_name submitter gbrs_paired_asm paired_asm_comp ftp_path excluded_from_refseq relation_to_type_material'.split()
genbank_summary                     = pd.read_table('/work/assembly_summary_genbank.txt', comment='#', header=None, names=header, dtype={'taxid':str, 'infraspecific_name':str})
genbank_summary['refseq_category']  = genbank_summary['refseq_category'].str.lower()
genbank_summary['assembly_level']   = genbank_summary['assembly_level'].str.lower()
genbank_summary['genome_rep']       = genbank_summary['genome_rep'].str.lower()
genbank_summary.set_index( 'assembly_accession', inplace=True )

refseq_summary                     = pd.read_table( '/work/assembly_summary_refseq.txt', comment='#', header=None, names=header, dtype={'taxid':str, 'infraspecific_name':str} )
refseq_summary['refseq_category']  = refseq_summary['refseq_category'].str.lower()
refseq_summary['assembly_level']   = refseq_summary['assembly_level'].str.lower()
refseq_summary['genome_rep']       = refseq_summary['genome_rep'].str.lower()
refseq_summary.set_index( 'assembly_accession', inplace=True )

assembly_summary = refseq_summary.append(genbank_summary)
assembly_summary.index = [index.replace('_', '').split('.')[0] for index in assembly_summary.index]

print 'yeah'
if os.path.isdir(sys.argv[1]):
    print 'Input is a folder'

    tree_folder = sys.argv[1]
    for tree_file in os.listdir(tree_folder):
        if not tree_file.endswith('.treefile') and not tree_file.endswith('.tree'):
            continue
        print tree_file

        tree = ete3.Tree('%s/%s' %(tree_folder, tree_file), format=1)
        out  = open('%s/%s.figTree' %(tree_folder, tree_file), 'wb')
        out.write("#NEXUS\nbegin taxa;\n\tdimensions ntax=%i;\n\ttaxlabels\n" %len(tree))
        for node in tree.traverse():
            if node.is_leaf():
                if node.name.startswith('GCF_') or node.name.startswith('GCA_'):
                    genome, gene = re.search('^(GC[FA]_\d+\.\d+)\|(\S+)$', node.name).groups()
                else:
                    genome, gene = re.search('^(\S+)_(\S+)$', node.name).groups()

                if genome not in assembly_summary.index:
                    out.write('\t%s\n' %(node.name))
                    continue

                out.write('\t%s ' %(node.name))
                comment = ['source_name="%s"' %assembly_summary.loc[genome, 'organism_name']]
                for rank in ['class', 'phylum', 'order', 'family']:
                    if rank in lineage:
                        comment.append('tax_%s="%s"' %(rank, assembly_summary.loc[genome, rank]))
                out.write('[&%s]\n' %' '.join(comment))

            else:
                if node.name:
                    if '/' in node.name:
                        write_format = 1
                        aLRT, UFBoot = node.name.split('/')
                        node.name = '[&UFBoot=%.2f,aLRT=%.2f]' %(float(UFBoot), float(aLRT))
                    else:
                        write_format = 0
                        node.support = float(node.name)

        newick_text = tree.write(format=write_format)
        if write_format:
            newick_text = re.sub('_&UFBoot_(\d+\.\d\d)_aLRT_(\d+\.\d\d)_', '[&UFBoot=\\1,aLRT=\\2]', newick_text)
        out.write(';\nend;\n')
        out.write('begin trees;\n\ttree tree_1 = [&R] %s\nend;' %newick_text)
        out.close()

else:
    print 'Input is a tree file'
    #
    # single tree
    print 'yeah'
    tree_folder = '.'
    tree_file   = sys.argv[1]
    tree = ete3.Tree('%s/%s' %(tree_folder, tree_file), format=1)
    out  = open('%s/%s.figTree' %(tree_folder, tree_file), 'wb')
    out.write("#NEXUS\nbegin taxa;\n\tdimensions ntax=%i;\n\ttaxlabels\n" %len(tree))
    branch_names = {}
    for node in tree.traverse():
        print 'yeah'
        if node.is_leaf():
            genome = node.name.split('_')[0]
            if genome == node.name:
                gene = ''
            else:
                gene = node.name.split('_')[1]
            if genome not in assembly_summary.index:
                out.write('\t%s\n' %(node.name))
                continue

            print genome

            taxid = assembly_summary.loc[genome, 'taxid']
            lineage = {j: i for i, j in ncbi.get_rank(ncbi.get_lineage(taxid)).items()}
            lineage_names = ncbi.get_taxid_translator(lineage.values())

            out.write('\t%s ' % (node.name))
            comment = ['source_name="%s"' %assembly_summary.loc[genome, 'organism_name']]
            for rank in ['class', 'phylum', 'order', 'family']:
                if rank in lineage:
                    comment.append('tax_%s="%s"' % (rank, lineage_names[lineage[rank]]))
            out.write('[&%s]\n' %' '.join(comment))

        else:
            if node.name:
                if '/' in node.name:
                    write_format = 1
                    aLRT, UFBoot = node.name.split('/')
                    node.name = '[&UFBoot=%.2f,aLRT=%.2f]' %(float(UFBoot), float(aLRT))
                else:
                    write_format = 0
                    node.support = float(node.name)

    newick_text = tree.write(format=write_format)
    if write_format:
        newick_text = re.sub('_&UFBoot_(\d+\.\d\d)_aLRT_(\d+\.\d\d)_', '[&UFBoot=\\1,aLRT=\\2]', newick_text)
    out.write(';\nend;\n')
    out.write('begin trees;\n\ttree tree_1 = [&R] %s\nend;' %newick_text)
    out.close()
