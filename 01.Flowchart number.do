*** 01. flowchart number
cd "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\data\"

*** read study_data.sas7bdat from S: .
import sas using "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\data\study_data.sas7bdat", clear
rename *, lower
compress
gen age=year(indexdate)-year(birthdate)
save study_data, replace

*** If there is new vars to be included, try to run this part and merge new vars
use controloutcomes, clear
compress
save, replace

merge 1:1 lopnr using study_data
drop gisymp* c19*
save study_data, replace


use study_data, clear
* COunts for flowchar
quietly count
local total=r(N)

quietly count if t2d_indication_ndr==1
local t2d=r(N)

quietly count if t2d_indication_ndr==1 & pre3y_3drug==1
local pre3ydrug=r(N)

quietly count if t2d_indication_ndr==1 & pre3y_3drug==0 & pre1y_pancreatitis==1
local pre1yap=r(N)

quietly count if t2d_indication_ndr==1 & pre3y_3drug==0 & pre1y_pancreatitis==0 & preindex_biliarycancer==1
local prebilcan=r(N)

quietly count if t2d_indication_ndr==1 & pre3y_3drug==0 & pre1y_pancreatitis==0 & preindex_biliarycancer==0 & preindex_pancreacancer==1
local prepanccan=r(N)

quietly count if t2d_indication_ndr==1 & pre3y_3drug==0 & pre1y_pancreatitis==0 & preindex_biliarycancer==0 & preindex_pancreacancer==0 & healthcare_region==""
local noregion=r(N)

quietly count if t2d_indication_ndr==1 & pre3y_3drug==0 & pre1y_pancreatitis==0 & preindex_biliarycancer==0 & preindex_pancreacancer==0 & healthcare_region!=""
local studypop=r(N)

quietly count if t2d_indication_ndr==1 & pre3y_3drug==0 & pre1y_pancreatitis==0 & preindex_biliarycancer==0 & preindex_pancreacancer==0 & healthcare_region!="" & expo=="glp1"
local glp1=r(N)
quietly count if t2d_indication_ndr==1 & pre3y_3drug==0 & pre1y_pancreatitis==0 & preindex_biliarycancer==0 & preindex_pancreacancer==0 & healthcare_region!="" & expo=="dpp4"
local dpp4=r(N)
quietly count if t2d_indication_ndr==1 & pre3y_3drug==0 & pre1y_pancreatitis==0 & preindex_biliarycancer==0 & preindex_pancreacancer==0 & healthcare_region!="" & expo=="sglt2"
local sglt2=r(N)


di "The source population comprised `total' individuals with at least one dispensation of any of the study drugs, of which `t2d' had a prior record of type 2 diabetes. From these, we excluded `pre3ydrug' with a prior dispensation of the study drugs, `pre1yap' previously hospitalised for acute pancreatitis, `prebilcan' with a previous diagnosis of biliary cancer, and `prepanccan' with a previous diagnosis of pancreatic cancer. We additionally excluded `noregion' with no information of the region they inhabit in. This yielded a total study population of `studypop' individuals (`=round(100*`studypop'/`total',1)'% retained). Of these, `glp1' were initiators of GLP-1 RAs, `dpp4' were initiators of DPP-4 inhibitors, and `sglt2' were initiators of SGLT-2 inhibitors (Figure 1). " 
