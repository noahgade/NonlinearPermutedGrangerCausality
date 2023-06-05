import numpy as np
import rdata
import cmlp

dat = rdata.conversion.convert(rdata.parser.parse_file("simdata1000b.RData"))

output = np.zeros((400, 7, 7), dtype = float)

for dat_indx in range(400):
    data = np.asmatrix(dat['simdata1000b'][dat_indx])
    dimension = data.shape[-1]
    output[dat_indx, :, :] = cmlp.fit_cmlpH(data, dimension, 3, 100, 0.5, 10000)
    
np.save("/home/ndg5e/NPGC3/cmlpH1000b.npy", output)