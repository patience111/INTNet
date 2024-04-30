import Bio.SeqIO as sio
import tensorflow as tf
import numpy as np
from sklearn.preprocessing import LabelBinarizer
from tensorflow.keras.utils import to_categorical
import random
import os
os.environ['CUDA_VISIBLE_DEVICES'] = '2'
import tqdm
#strategy = tf.distribute.MirroredStrategy()

#load model
filterm = tf.keras.models.load_model(os.path.join(os.path.dirname(__file__), '../model/inti-aels_tall.h5'))
classifier = tf.keras.models.load_model(os.path.join(os.path.dirname(__file__), '../model/inti-intitypes-ls_tall.h5'))
bac_host = tf.keras.models.load_model(os.path.join(os.path.dirname(__file__), '../model/inti-host-ls_tall.h5'))
asso_args = tf.keras.models.load_model(os.path.join(os.path.dirname(__file__), '../model/inti-argls_tall.h5'))

#encode, encode all the sequence to 1600 aa length
char_dict = {}
chars = 'ACDEFGHIKLMNPQRSTVWXYBJZ'
new_chars = "ACDEFGHIKLMNPQRSTVWXY"
for char in chars:
    temp = np.zeros(22)
    if char == 'B':
        for ch in 'DN':
            temp[new_chars.index(ch)] = 0.5
    elif char == 'J':
        for ch in 'IL':
            temp[new_chars.index(ch)] = 0.5
    elif char == 'Z':
        for ch in 'EQ':
            temp[new_chars.index(ch)] = 0.5
    else:
        temp[new_chars.index(char)] = 1
    char_dict[char] = temp

def encode(seq):
    char = 'ACDEFGHIKLMNPQRSTVWXY'
    train_array = np.zeros((1600, 22))
    for i in range(1600):
        if i<len(seq):
            train_array[i] = char_dict[seq[i]]
        else:
            train_array[i][21] = 1
    return train_array

def test_encode(tests):
    tests_seq = []
    for test in tests:
        tests_seq.append(encode(test))
    tests_seq = np.array(tests_seq)
    
    return tests_seq

def newEncodeVaryLength(seq):
    char = 'ACDEFGHIKLMNPQRSTVWXY'
    mol = len(seq) % 16
    dimension1 = len(seq) - mol + 16
    train_array = np.zeros((dimension1,22))
    for i in range(dimension1):
        if i < len(seq):
            train_array[i] = char_dict[seq[i]]
        else:
            train_array[i][21] = 1
    
    return train_array

def test_newEncodeVaryLength(tests):
    tests_seq = []
    for test in tests:
        tests_seq.append(newEncodeVaryLength(test))
    tests_seq = np.array(tests_seq)
    
    return tests_seq

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def filter_prediction_batch(seqs):
    predictions = []
   # for seq in seqs:
    #    temp = model.predict(np.array([seq]))
     #   predictions.append(temp)
    temp = filterm.predict(seqs, batch_size = 512)
    predictions.append(temp)
    return predictions

def prediction(seqs):
    predictions = []
    for seq in seqs:
        temp = model.predict(np.array([seq]))
        predictions.append(temp)
    return predictions

def reconstruction_simi(pres, ori):
    simis = []
    reconstructs = []
    argmax_pre = np.argmax(pres[0], axis=2)
    for index, ele in enumerate(argmax_pre):
        length = len(ori[index])
        count_simi = 0
        #reconstruct = ''
        if length >= 1600:
            align = 1600
        else:
            align = length
        for pos in range(align):
            if chars[ele[pos]] == ori[index][pos]:
                count_simi += 1
            #reconstruct += chars[np.argmax(ele[pos])]
        simis.append(count_simi / length)
        #reconstructs.append(reconstruct)
    return simis

intis_labels = ['Group09', 'Group10', 'Group13', 'Group03', 'Group02', 'Group01', 'Group05', 'Group11',
 'Group06', 'Unclassified', 'Group12', 'Group04', 'Group07', 'Group14', 'Group15', 'Group08']

hosts_labels = ['Gammaproteobacteria',  'Unassigned',  'Betaproteobacteria',  'Deltaproteobacteria',
           'Planctomycetia',  'Spirochaetia',  'Anaerolineae',  'Verrucomicrobiae',  'Opitutae',
           'Epsilonproteobacteria',  'Ignavibacteria',  'Gemmatimonadetes',  'Candidatus Brocadiia',
           'Balneolia',  'Spartobacteria',  'Kiritimatiellae',  'Blastocatellia',  'Chlorobia',  'Nitrospira',
           'Actinomycetia',  'Hydrogenophilalia',  'Phycisphaerae',  'Chloroflexia',  'Alphaproteobacteria',
           'Holophagae',  'Acidithiobacillia',  'Zetaproteobacteria',  'Coriobacteriia',  'Lentisphaeria',
           'Rhodothermia',  'Caldilineae',  'Bacteroidia',  'Acidobacteriia',  'Clostridia',  'Calditrichia',
           'Thermoanaerobaculia',  'Thermodesulfobacteria',  'Vicinamibacteria',  'Deinococci',  'Bacilli',
           'Chitinophagia',  'Ardenticatenia',  'Thermoflexia',  'Acidimicrobiia',  'Flavobacteriia',
           'Tichowtungiia',  'Thermodesulfovibrionia',  'Fimbriimonadia',  'Chitinivibrionia',  'Oligoflexia',
           'Candidatus Fermentibacteria (class)',  'Candidatus Lambdaproteobacteria',
           'Candidatus Muproteobacteria',  'Abditibacteria',  'Candidatus Binatia',  'Limnochordia',
           'Synergistia',  'Candidatus Thermofonsia',  'Thermoleophilia',  'Cytophagia',  'Tepidiformia',
           'Aquificae',  'Chrysiogenetes',  'Gloeobacteria',  'Fibrobacteria',  'Chitinispirillia',
           'Dehalococcoidia',  'Candidatus Ozemobacteria',  'Candidatus Polarisedimenticolia',
           'Erysipelotrichia']

args_labels = ['MLS', 'aminocoumarin', 'aminoglycoside', 'bacitracin', 'beta-lactam', 'bleomycin',
        'chloramphenicol', 'elfamycin', 'ethambutol', 'fosfomycin', 'fosmidomycin', 'fusidic_acid',
        'glycopeptide', 'isoniazid', 'kasugamycin', 'multidrug', 'mupirocin', 'nitrofurantoin',
        'nitroimidazole', 'peptide', 'pleuromutilin', 'polymyxin', 'pyrazinamide', 'qa_compound',
        'quinolone', 'rifamycin', 'streptothricin', 'sulfonamide', 'tetracenomycin', 'tetracycline',
        'triclosan', 'trimethoprim', 'tunicamycin']

intis_prepare = sorted(intis_labels)
intis_dic = {}
for index, ele in enumerate(intis_prepare):
    intis_dic[index] = ele


hosts_prepare = sorted(hosts_labels)
hosts_dic = {}
for index, ele in enumerate(hosts_prepare):
    hosts_dic[index] = ele


args_prepare = sorted(args_labels)
#with strategy.scope():
def intnet_lsaa(input_file, outfile):
    #cuts: 0.4474056985591711,  0.44070125025206697, 0.46015594589334374, 0.482238570677877, 0.4825442827997788    
    cut = 0.46260914963644756
    print('reading in test file...')
    test = [i for i in sio.parse(input_file, 'fasta')]

    with open(os.path.join(os.path.dirname(__file__), "../results/" + outfile) , 'a') as f:
            f.write('test_id' + '\t' + 'inti_type' + '\t' +  'pre_prob' + '\t' + 'bacterial_host' + '\t' + \
                    'pre_prob' + '\t' + 'resistance_category' + '\n')

    print('encoding test file...')
    print('reconstruct, simi...')
    for idx, test_chunk in enumerate(list(chunks(test, 10000))):
    #test_ids = [ele.id for ele in test]
        testencode = test_encode(test_chunk)
        testencode_pre = filter_prediction_batch(testencode) # if huge volumn of seqs (~ millions) this will be change to create batch in advance 
        simis = reconstruction_simi(testencode_pre, test_chunk)
    #results = calErrorRate(simis, cut) 
    #passed = []
        passed_encode = [] ### notice list and np.array
        passed_idx = []
        notpass_idx = []
        for index, ele in enumerate(simis):
            if ele >= cut:
                #passed.append(test[index])
                passed_encode.append(testencode[index])
                passed_idx.append(index)
            else:
                notpass_idx.append(index)
    
        ###classification
        print('classifying...')
        
        if len(passed_encode) > 0:
            classifications = classifier.predict(np.stack(passed_encode, axis=0), batch_size = 512)
            classification_argmax = np.argmax(classifications, axis=1)
            classification_max = np.max(classifications, axis=1)
            
            hosts = bac_host.predict(np.stack(passed_encode, axis=0), batch_size = 512)
            hosts_argmax = np.argmax(hosts, axis=1)
            hosts_max = np.max(hosts, axis=1)
            
            args = asso_args.predict(np.stack(passed_encode, axis=0), batch_size = 512)
            args1 = np.squeeze(np.where(args >= 0.5, 1, 0))
            #args2 = np.array(args_prepare)[args1]
            print('args1.shape: ', args1.shape)
            
            if len(args1.shape) == 1:
                args1 = args1.reshape((len(args1.shape), len(args_prepare)))
            
            args2 = [[x for x, pred in zip(args_prepare, y) if pred ==1] for y in args1]
#            args2 = []
#            for y in args1:
#                print('type y :', type(y))
#                print('len y: ', len(y))
#                for x, pred in zip(args_prepare, y):
#                    if pred == 1:
#                        args2.append(x)

            assert len(args2)==len(passed_idx)
            
            inti = {}
            host = {}
            args2_e1 = {}
            
            for i, ele in enumerate(passed_idx):
                inti[ele] = [classification_max[i], intis_dic[classification_argmax[i]]]
                host[ele] = [hosts_max[i], hosts_dic[hosts_argmax[i]]]
                args2_e1[ele] = args2[i]

            ### output
            print('writing output...')
            with open(os.path.join(os.path.dirname(__file__), "../results/" + outfile) , 'a') as f:
                for idx, ele in enumerate(test_chunk):
                    if idx in passed_idx:
                        f.write(test_chunk[idx].id + '\t')
                        f.write(inti[idx][-1] + '\t')
                        f.write(str(inti[idx][0]) + '\t')
                        f.write(host[idx][-1] + '\t')
                        f.write(str(host[idx][0]) + '\t')
                        f.write(','.join(args2_e1[idx]) + '\n')

                    if idx in notpass_idx:
                        f.write(test_chunk[idx].id + '\t')
                        f.write('non-inti' + '\t' + '' + '\t' + '' + '\t' + '' + '\t' + '' + '\n')
            
        if len(passed_encode) == 0:
            print('no seq passed!')
            pass

