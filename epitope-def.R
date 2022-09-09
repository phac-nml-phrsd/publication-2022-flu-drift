###
###   DEFINITION OF EPITOPES POSITIONS 
###


#' Position of H3N2 epitopes according to several definitions
#' 
get_epitope_positions_H3N2 <- function() {
  
  sites.narrow = c(60:64, 66, 67, 69, 70, 73, 75, 78, 79, 83, 91, 94, 96:99, 
                   102:104, 107, 108, 110, 112, 118, 119, 125, 133, 137, 138, 
                   140, 142, 144:149, 151, 153, 154, 156, 158:162, 166, 168, 
                   171:176, 179:181, 183, 184, 186:193, 195, 198, 
                   202:206, 208:210, 212:214, 217, 219, 223:225, 228:235, 
                   242:246, 254, 256, 258, 260, 262:264, 276:278, 281, 
                   289, 291, 292, 294:296, 
                   310, 313, 315, 316, 320, 321, 323:328 )
  
  # approximately +1AA buffer from `narrow`
  sites.broad = c(59:80, 82:84, 90:113, 117:120, 124:126, 132:134, 136:163, 
                  165:199, 201:220, 222:236, 241:247, 253:265, 275:282, 
                  288:297, 309:317, 319:329 )
  
  # approximately +3AA buffer from `narrow`
  sites.broad3 = c(57:86, 88:128, 130:249, 251:267, 273:284, 286:299, 307:331)
  
  # From Smith et al. 2010
  sites.smith2010 = c(41, 66, 69, 70, 78, 91, 98, 99, 138, 140, 147, 149, 153, 
                      159:162, 171, 172, 174, 176, 180, 188, 190, 204:206, 209, 
                      212, 213, 217, 218, 223, 229, 233, 238, 241, 246, 260, 
                      276, 278, 292, 294)
  
  # From Koele et al. 2013 and Beer et al. 2018
  
  sites.koele.beer = c(161, 171, 172, 174, 175, 205, 206, 209)
  
  sites.full.HA1  = c(17:345)
  
  sites.full.HA12 = c(17:566)
  
  
  sites.list = list(
    narrow     = sites.narrow,
    broad      = sites.broad,
    broad3     = sites.broad3,
    misc1      = sites.smith2010,
    misc2      = sites.koele.beer,
    full.HA1   = sites.full.HA1,
    full.HA12  = sites.full.HA12
  )
  return(sites.list)
}


get_epitope_positions_H1N1 <- function() {
  
  sites.narrow = c(87:92, 141, 142, 154:159, 170:174, 176:181, 183:187, 
                   201:212, 220:222, 238, 239, 252:254) 
  
  # approximately +1AA buffer from `narrow`
  sites.broad = c(86:93, 140:143, 153:160, 169:188, 200:213, 
                  219:223, 237:240, 251:255) 
  
  # approximately +3AA buffer from `narrow`
  sites.broad3 = c(84:95, 138:145, 151:162, 167:190, 198:215, 
                   217:225, 235:242, 249:257)
  
  # Combined - Koel (2013), Koel (2015), Rudneva (2012), 
  # Manicassamy (2010), O'Donnell (2012)
  sites.combined = c(136, 142, 144, 147, 148, 158, 169:173, 180, 200, 204, 241)
  
  sites.full.HA1  = c(18:344)
  
  sites.full.HA12 = c(18:566)
  
  sites.list = list(
    narrow    = sites.narrow,
    broad     = sites.broad,
    broad3    = sites.broad3,
    misc1     = sites.combined,
    misc2     = sites.combined,
    full.HA1  = sites.full.HA1,
    full.HA12 = sites.full.HA12
  )
  return(sites.list)
}

get_epitope_positions_B <- function() {
  
  sites.narrow = c(63, 71, 88:94, 131:152, 156:165, 177:185,
                   197:199, 212:220, 244:259) 
  
  sites.broad = c(62:64, 70:72, 87:95, 130:153, 155:166, 
                  176:186, 196:200, 211:221, 243:260) 
  
  # approximately +3AA buffer from `narrow`
  sites.broad3 = c(60:66, 68:74, 85:97, 128:168, 174:188, 
                   194:202, 209:223, 241:262)
  
  # Site Koel et al 2013
  sites.koele2013 = c(180,181)
  
  sites.full.HA1  = c(16:360)
  
  sites.full.HA12 = c(16:585)
  
  
  sites.list = list(
    narrow    = sites.narrow,
    broad     = sites.broad,
    broad3    = sites.broad3,
    misc1     = sites.koele2013,
    misc2     = sites.koele2013,
    full.HA1  = sites.full.HA1,
    full.HA12 = sites.full.HA12
  )
  return(sites.list)
}



get_epitope_positions <- function(virus) {
  
  if(virus=='H1N1') return(get_epitope_positions_H1N1())
  if(virus=='H3N2') return(get_epitope_positions_H3N2())
  if(virus=='B')    return(get_epitope_positions_B())
  
}


