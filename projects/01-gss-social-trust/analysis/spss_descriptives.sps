* GSS 2024 descriptive research syntax. Prepared but not run in IBM SPSS.
* Use actual survey design software for survey-correct standard errors.
GET FILE='data/raw/GSS2024.sav'.
DATASET NAME gss2024 WINDOW=FRONT.
SELECT IF (SAMPLE=13).
EXECUTE.
DISPLAY DICTIONARY.
FREQUENCIES VARIABLES=TRUST HAPPY CONFINAN CONEDUC WTSSNRPS WTSSPS
 /ORDER=ANALYSIS.
RECODE TRUST (1=1) (2=2) (3=3) (ELSE=SYSMIS) INTO trust_valid.
RECODE HAPPY (1=1) (2=2) (3=3) (ELSE=SYSMIS) INTO happy_valid.
RECODE CONFINAN (1=1) (2=2) (3=3) (ELSE=SYSMIS) INTO confinan_valid.
RECODE CONEDUC (1=1) (2=2) (3=3) (ELSE=SYSMIS) INTO coneduc_valid.
VALUE LABELS trust_valid 1 'Most people can be trusted' 2 'Cannot be too careful' 3 'Depends'.
VALUE LABELS happy_valid 1 'Very happy' 2 'Pretty happy' 3 'Not too happy'.
VALUE LABELS confinan_valid coneduc_valid 1 'Great deal' 2 'Only some' 3 'Hardly any'.
COMPUTE trust_seen = NOT MISSING(trust_valid).
COMPUTE confinan_seen = NOT MISSING(confinan_valid).
COMPUTE coneduc_seen = NOT MISSING(coneduc_valid).
COMPUTE happy_seen = NOT MISSING(happy_valid).
EXECUTE.
CROSSTABS /TABLES=trust_seen BY confinan_seen /CELLS=COUNT ROW COLUMN.
CROSSTABS /TABLES=trust_seen BY coneduc_seen /CELLS=COUNT ROW COLUMN.
CROSSTABS /TABLES=trust_seen BY happy_seen /CELLS=COUNT ROW COLUMN.
FREQUENCIES VARIABLES=trust_valid happy_valid confinan_valid coneduc_valid.
CROSSTABS /TABLES=trust_valid BY happy_valid /CELLS=COUNT ROW.
COMPUTE eligibleweight = WTSSNRPS.
IF (eligibleweight <= 0) eligibleweight=$SYSMIS.
WEIGHT BY eligibleweight.
FREQUENCIES VARIABLES=trust_valid happy_valid confinan_valid coneduc_valid.
CROSSTABS /TABLES=trust_valid BY happy_valid /CELLS=COUNT ROW.
WEIGHT OFF.
* This syntax will generate genuine SPSS Viewer results only when executed in SPSS.
