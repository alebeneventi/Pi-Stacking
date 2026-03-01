
dscf = 0
elst = 0
exch = 0
ind = 0
exch_ind = 0
disp = 0
sapt_total = 0
xsapt_total = 0

#``dhf_name`` is the name of the output file from the delta SCF job
dhf_name = 'dhf.in.out'
with open(dhf_name, 'r') as dhf:
    for line in dhf:
        val = 'Delta(SCF)'
        if line.find(val) !=-1 and dscf == 0:
            #Find the first instance of ''Delta(SCF)'' from the output file.
            #This value is in kcal/mol
            dscf = float(line.split()[1])

#``sapt_name`` is the name of the output file from the SAPT job without charge embedding job
sapt_name = 'sapt.in.out'
with open(sapt_name, 'r') as dhf:
    for line in dhf:
        elst_str = 'E1_elst'
        exch_str = 'E1_exch'
        ind_str = 'E2_ind'
        exch_ind_str = 'E2_exch-ind'
        disp_str = 'MBD disperion'
        sapt_total_str = 'SAPT corrected total energy'
        if line.find(elst_str) !=-1 and elst == 0:
            #Find the first instance of ''E1_elst'' from the output file.
            #This value is in kcal/mol
            elst = float(line.split()[1])
        if line.find(exch_str) !=-1 and exch == 0:
            #Find the first instance of ''E1_exch'' from the output file.
            #This value is in kcal/mol
            exch = float(line.split()[1])
        if line.find(ind_str) !=-1 and ind == 0:
            #Find the first instance of ''E2_ind'' from the output file.
            #This value is in kcal/mol
            ind = float(line.split()[1])
        if line.find(exch_ind_str) !=-1 and exch_ind == 0:
            #Find the first instance of ''E2_exch'' from the output file.
            #This value is in kcal/mol
            exch_ind = float(line.split()[1])
        if line.find(disp_str) !=-1 and disp == 0:
            #Find the first instance of ''MBD disperion'' from the output file.
            #This value is in kcal/mol
            disp = float(line.split()[2])
        if line.find(sapt_total_str) !=-1 and sapt_total == 0:
            #Find the first instance of ''SAPT corrected total energy'' from the output file.
            #This value is in Hartrees
            sapt_total = float(line.split()[5])

#``sapt_em_name`` is the name of the output file from the SAPT job with charge embedding job
sapt_em_name = 'xsapt.in.out'
with open(sapt_em_name, 'r') as dhf:
    for line in dhf:
        val = 'SAPT corrected total energy'
        if line.find(val) !=-1 and xsapt_total == 0:
            #Find the first instance of ''SAPT corrected total energy'' from the output file.
            #This value is in Hartrees
            xsapt_total = float(line.split()[5])

total = elst+exch+ind+exch_ind+(xsapt_total-sapt_total)*627.5+disp+dscf
print("Elst,Exch,Ind,Disp,Total")
print('{},{},{},{},{}'.format(elst,exch,round(ind+exch_ind+dscf+(xsapt_total-sapt_total)*627.5, 3),disp,round(total, 3)))
#elst is the first-order electrostatic term from the SAPT job without charge embedding
#exch is the first-order exchange term from the SAPT job without charge embedding
#ind is sum of the second-order induction and exchange-induction from the SAPT job without charge embedding,
#the delta SCF value from the dhf job, and the difference in the SAPT corrected total energy between the 
#SAPT job with and without charge embedding.
#disp is the MBD dispersion value