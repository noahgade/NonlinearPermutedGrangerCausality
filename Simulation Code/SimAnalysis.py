import numpy as np
import pandas as pd

GL250 = cmlpGL250[:, 0, :]
pd.DataFrame(GL250).to_csv("GL250.csv")

GL500 = cmlpGL500[:, 0, :]
pd.DataFrame(GL500).to_csv("GL500.csv")

GL1000 = cmlpGL1000[:, 0, :]
pd.DataFrame(GL1000).to_csv("GL1000.csv")


GSGL250 = cmlpGSGL250[:, 0, :]
pd.DataFrame(GSGL250).to_csv("GSGL250.csv")

GSGL500 = cmlpGSGL500[:, 0, :]
pd.DataFrame(GSGL500).to_csv("GSGL500.csv")

GSGL1000 = cmlpGSGL1000[:, 0, :]
pd.DataFrame(GSGL1000).to_csv("GSGL1000.csv")


H250 = cmlpH250[:, 0, :]
pd.DataFrame(H250).to_csv("H250.csv")

H500 = cmlpH500[:, 0, :]
pd.DataFrame(H500).to_csv("H500.csv")

H1000 = cmlpH1000[:, 0, :]
pd.DataFrame(H1000).to_csv("H1000.csv")


GL1000b = cmlpGL1000b[:, 0, :]
pd.DataFrame(GL1000b).to_csv("GL1000b.csv")

GSGL1000b = cmlpGSGL1000b[:, 0, :]
pd.DataFrame(GSGL1000b).to_csv("GSGL1000b.csv")

H1000b = cmlpH1000b[:, 0, :]
pd.DataFrame(H1000b).to_csv("H1000b.csv")