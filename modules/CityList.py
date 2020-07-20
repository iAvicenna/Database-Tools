"""
@author: Sina Tureli

A module for organizing cities into classes. Some functionalities
are for cities like HONG_KONG, it generates aliases like
HONGKONG, HONG-KONG, HONG KONG. This helps match antigen and
sera ids more precisely. Moreover when trying to match a
query like HONG KONG to a city named HONGKONG, aliasing
is also applied to the query (not just to the city record itself). The search
supports regular expressions.

Searching functionality also allows searching for a city by edit distance
so searches for city names with typos like HONK KONG can still find the correct
city. 

"""

import editdistance as ed
import numpy as np
import re

def generate_city_aliases(name):   
    aliases = []
    aliases.append(name.upper())
    
    replace_chars = ['-', '_', ' ', '/']
    
    for ind1, char1 in enumerate(replace_chars):
        for ind2, char2 in enumerate(replace_chars[ind1+1:]+[''], start=ind1+1):
            
            if char1 in aliases[0]:
                aliases.append(aliases[0].replace(char1, char2))

    return aliases

class City:
    
    def __init__(self, name, abb):
        
        self.name = name.upper()
        self.abb = abb
        self.aliases = generate_city_aliases(name)
        
          
    def match(self, query, aliasing=True, exact_match=False):

        '''Try to match query to city name
       
        If aliasing = True, then aliases for query and city name is
        taken onto account
        '''
        
        assert isinstance(query, str), 'query must be a string'
        
        if aliasing:
            query_aliases =  generate_city_aliases(query)
            city_aliases = self.aliases
        else: 
            query_aliases = [query]
            city_aliases = [self.name]
        
        for query in query_aliases: 
            if exact_match:
                query = '^' + query + '$'
            regexp = re.compile(query, re.IGNORECASE)
            if any([regexp.search(x) for x in city_aliases]):
                return True
          
        return False
        
        
    def ematch(self, query, aliasing=True):
        
        '''
        Match query to aliases via editdistance
        '''
        
        if aliasing:
            lookup_vals = self.aliases
        else:
            lookup_vals = [self.name]
        
        edistances = [ed.eval(x, query) for x in lookup_vals]
        return min(edistances)
                    
    
    def __str__(self):
        
        return 'City {} with abbreviation {}'.format(self.name, self.abb)
    
    def __repr__(self):
        
        return '{}:{}'.format(self.abb, self.name)
    
class CityList:
    
    ''' List of Cities:
    
    It is initialized from the city_map dictionary at the end of this module.
    
    It counts nonunique cities, abbreviations and problematic
    (abbreviation, city name) pairs. It also does a health
    check which reports non uniqueness problems etc.
    
    Contains search functionalities to see whether if a given query
    exists in this list or not. The search functionality is built 
    from the search functionalities of the City class above itself
    so allows aliasing and edit distance searches.
    
    '''
    
    def __init__(self):
        
        self.cities = []
        
        self.nonunique_cities = []
        self.nonunique_abbs = []
        self.problematic_names = []
        
        for abb in city_map:
            city_name = city_map[abb]
                
            if len(abb) ==0 or len(abb) > 2:
               print('Abbreviation {} for city {} too short or too long.'.format(abb, city_name))
               self.problematic_names.append((abb,city_name)) 
            if len(city_name) < 3:
               print('Warning: City name {} with abbreviation {} is too short.'.format(city_name, abb))
               self.problematic_names.append((abb,city_name))
            
            self.add_city(City(city_name,abb))
    
        self.health_check()
    
    def add_city(self, city):
        
        assert isinstance(city, City)
        
        self.cities.append(city)

     
    def abb_to_city(self, abb):
        
        I = [ind for ind,city in enumerate(self.cities) if city.abb == abb]
        
        return [self.cities[x] for x in I]
    
    def name_to_city(self, name):
        
        I = [ind for ind,city in enumerate(self.cities) if city.name == name]
        
        return [self.cities[x] for x in I]
        
    def search(self, query, split=False, aliasing=False, exact_match=False):
        
        assert isinstance(query, str), 'Query to search should be a string'
        
        if split:
            query_list = query.split('/')
            query_list = [x for x in query_list if len(x)>2]
        else:
            query_list = [query]
            
        
        I = []
        
        for q in query_list:
            I += [ind  for ind,city in enumerate(self.cities) if city.match(q, aliasing=aliasing, exact_match=exact_match)]

        
        return [self.cities[ind] for ind in I]
    
    
    def esearch(self, query, aliasing=True):
        
        '''
        Search via edit distance. Return the matches with minimum edit distance
        if this edit distance is lower than half the length of the query.
        '''
        
        assert isinstance(query, str), 'Query to esearch should be a string'
        
        edistances = np.array([city.ematch(query, aliasing=aliasing) for city in self.cities])
        
        min_edist = np.min(edistances)
        
        if min_edist < len(query)/2:
            I = np.argwhere(edistances == min_edist).flatten()
            
            return ([self.cities[x] for x in I], min_edist)
        else:
            
            return ([], len(query)/2)
        
    def health_check(self):
        '''
        Check for non-unique city names and abbreviations.
        '''
        
        self.nonunique_cities = []
        for city in self.cities:
            
            I = self.name_to_city(city.name)
            
            if len(I)>1 and city.name not in self.nonunique_cities:
                self.nonunique_cities.append(city.name)
                print('City {} appears more than once:'.format(city.name))
                print(I)
        
        
        abbs = [x.abb for x in self.cities]
        self.nonunique_abbs = []
        
        for abb in abbs:
            
            I = abbs.count(abb)
            
            if I>1 and abb not in self.nonunique_abbs:
                self.nonunique_abbs.append(abb)
                print('Abbreviation {} appears {} times'.format(abb, I))
                print(self.abb_to_city(abb))
        
    def __str__(self):
        
        return '''{} cities with {} non-unique cities, {} non-unique 
                abbreviations and {} problematic names'''.format(len(self.cities),
                len(self.nonunique_cities),len(self.nonunique_abbs),
                len(self.problematic_names))
   
    
city_map = {
    "-A": "SAUDI-ARABIA",
    "-B": "SUPHAN-BURI",
    "-D": "NORDRHEIN-WESTFALEN",
    "-G": "FRENCH-GUIANA",
    "-H": "SENDI-H",
    "-I": "RHODE-ISLAND",
    "-K": "SA-KAEO",
    "-L": "SRI-LANKA",
    "-M": "PUERTO-MONTT",
    "-N": "SAKON-NAKHON",
    "-O": "DISTRICT-OF-COLUMBIA",
    "-P": "ST-PETERSBURG",
    "-S": "SAN-SEBASTIAN",
    "-T": "SAN-ANTONIO",
    "-V": "MECKLENBURG-VORPO.",
    "AA": "ANNARBOR",
    "AB": "ALABAMA",
    "AC": "AICHI",
    "AD": "ALBANIA",
    "AE": "ALGERIA",
    "AF": "SOUTH-AFRICA",
    "AG": "ARGENTINA",
    "AH": "AUSTRIA",
    "AI": "ASTURIAS",
    "AJ": "SAKAI",
    "AK": "AKITA",
    "AL": "AUSTRALIA",
    "AM": "AMSTERDAM",
    "AN": "ANHUI",
    "AO": "AOMORI",
    "AP": "CHIANGMAI",
    "AQ": "ANG-THONG",
    "AR": "ARIZONA",
    "AS": "ALICESPRINGS",
    "AT": "ATLANTA",
    "AU": "AUCKLAND",
    "AV": "SANTIAGAO",
    "AX": "ALASKA",
    "AY": "AYTTHAYA",
    "AZ": "AYATTHAYA",
    "BA": "BANGKOK",
    "BB": "BILBAO",
    "BC": "BUCHAREST",
    "BD": "BADEN-WURTTEMBERG",
    "BE": "BEIJING",
    "BF": "BALEARES",
    "BG": "BELGIUM",
    "BH": "BELGRADE",
    "BI": "BILTHOVEN",
    "BJ": "BEIIJNG",
    "BK": "BANGKOKI",
    "BL": "BELEM",
    "BM": "BREMEN",
    "BN": "BERLIN",
    "BO": "BARCELONA",
    "BP": "BURIRUM",
    "BQ": "PATTANI",
    "BR": "BRISBANE",
    "BS": "BRASOV",
    "BT": "BAGKOK",
    "BU": "BUSAN",
    "BV": "PHARE",
    "BW": "PHICHIT",
    "BX": "BURIRAM",
    "BY": "BAYERN",
    "BZ": "BRAZIL",
    "CA": "CANBERRA",
    "CB": "CHIBA",
    "CC": "CHRISTCHURCH",
    "CD": "COLLINDALE",
    "CE": "CHEONNAM",
    "CF": "CALIFORNIA",
    "CG": "CHUNGNAM",
    "CH": "CHINA",
    "CI": "CAEN",
    "CJ": "CHEJU",
    "CK": "CHEONBUK",
    "CL": "CHILE",
    "CM": "CHIANGMAI",
    "CN": "CONNECTICUT",
    "CO": "COLORADO",
    "CP": "CLUJ",
    "CQ": "CANARIAS",
    "CR": "CHIANGRAI",
    "CS": "CALARASI",
    "CT": "CHITA",
    "CU": "CHANTABURI",
    "CV": "CARTAGENA",
    "CW": "CHANGWON",
    "CZ": "CZECHOSLOVAKIA",
    "DA": "DAEGU",
    "DB": "NY",
    "DC": "CA",
    "DD": "DUNDIN",
    "DE": "DELFT",
    "DF": "KYUNGGI",
    "DG": "PAU",
    "DH": "MILAN",
    "DI": "SALAMANCA",
    "DJ": "DAEJEON",
    "DK": "DAKAR",
    "DL": "CASABLANCA",
    "DM": "MEKNES",
    "DN": "DENMARK",
    "DO": "MARRAKECH",
    "DP": "AGADIR",
    "DQ": "ATHENS",
    "DR": "SANTANDER",
    "DS": "DESGENETTES",
    "DU": "DUNEDIN",
    "DV": "DEVA",
    "DW": "DARWIN",
    "DX": "TENNESEE",
    "DY": "BRANDYS",
    "DZ": "TRANG",
    "EA": "SIENA",
    "EB": "ANGTHONG",
    "EC": "ECUADOR",
    "ED": "NIEDERSACHSEN",
    "EE": "NEPAL",
    "EF": "KALASIN",
    "EG": "EGYPT",
    "EH": "EHIME",
    "EI": "EINDHOVEN",
    "EJ": "PHITSANULOK",
    "EK": "EKATERINBURG",
    "EL": "EL-SALVADOR",
    "EM": "PRACHUAPKHIRIKKAN",
    "EN": "ENGLAND",
    "EO": "CHACHOENGSAO",
    "EP": "NEW-HAMPSHIRE",
    "EQ": "EQUADOR",
    "ER": "TOMSK",
    "ES": "ENSCHEDE",
    "ET": "EXTREMAD",
    "EU": "EUSKADI",
    "EV": "PONTEVEORA",
    "EX": "EXTREMADURA",
    "EY": "KYOTO",
    "EZ": "EKATEZINBURG",
    "FA": "FATICK",
    "FB": "FUKUI",
    "FC": "FUKUOKA-C",
    "FD": "FES",
    "FE": "FIRENZE",
    "FF": "ANNECY",
    "FG": "MACU",
    "FH": "FOSHAN",
    "FI": "FINLAND",
    "FJ": "FIJI",
    "FK": "FUOKA",
    "FL": "FLORIDA",
    "FM": "SUKHBAATAR",
    "FN": "SAINSHAND",
    "FO": "FLORENCE",
    "FP": "CAMBODIA",
    "FQ": "HOLLAND",
    "FR": "FRANCE",
    "FS": "FUKUSHIMA",
    "FT": "DUBAI",
    "FU": "FUJIAN",
    "FV": "GANSUCHENGGUAN",
    "FW": "HONDURAS",
    "FX": "KALININGRAD",
    "FY": "ASTRAKHAN",
    "FZ": "FUZHOU",
    "GA": "GUAM",
    "GB": "GENOA",
    "GC": "GREECE",
    "GD": "GUANGDONG",
    "GE": "GENEVA",
    "GF": "GUILDFORD",
    "GG": "GEORGIA",
    "GH": "GIRONA",
    "GI": "GIFU",
    "GJ": "GUADALAJARA",
    "GK": "GOTEBURG",
    "GL": "GUALDALAGARA",
    "GM": "GUMMA",
    "GN": "GRANADA",
    "GO": "GOTENBORG",
    "GP": "GIFU-C",
    "GQ": "GALICIA",
    "GR": "GRONINGEN",
    "GS": "GANSU",
    "GT": "GOTEBORG",
    "GU": "GUIZHOU",
    "GV": "SEGOVIA",
    "GW": "GUNMA",
    "GX": "GUANGXI",
    "GY": "GERMANY",
    "GZ": "GUANGZHOU",
    "HA": "HAWAII",
    "HB": "HARBIN",
    "HC": "SACHSEN",
    "HD": "HOKKAIDO",
    "HE": "HEBEI",
    "HF": "HIROSHIMA",
    "HG": "HUNGARY",
    "HH": "SHIZUOKA",
    "HI": "HUBEI",
    "HJ": "HYOGO",
    "HK": "HONGKONG",
    "HL": "HIME",
    "HM": "HENAN",
    "HN": "HUNAN",
    "HO": "HOUSTON",
    "HP": "HONG-KONG",
    "HQ": "HAINAN",
    "HS": "HIROSHIMA-C",
    "HT": "HAMAMATU-C",
    "HU": "HUNTINGTON",
    "HV": "HANNOVER",
    "HW": "NAKHON-SAWAN",
    "HX": "SHAANXI",
    "HY": "CHANTHABURI",
    "HZ": "SHIZUOKA-C",
    "IA": "IASI",
    "IB": "IBARAKI",
    "IC": "INCHEON",
    "ID": "INDONESIA",
    "IE": "ICELAND",
    "IF": "IRAN",
    "IH": "IDAHO",
    "II": "INDIA",
    "IJ": "BEIJINGXUANWU",
    "IK": "ISHIKAWA",
    "IL": "ILLINOIS",
    "IM": "SHIMANE",
    "IN": "INDIANA",
    "IO": "SOLOMON-ISLANDS",
    "IP": "MURMANSK",
    "IQ": "STPETERSBURG",
    "IR": "IRELAND",
    "IS": "ISRAEL",
    "IT": "ITALY",
    "IU": "ISTANBUL",
    "IV": "INVERNESS",
    "IW": "IWATE",
    "IX": "BEIJINXUANWU",
    "IY": "CHAIYAPUM",
    "IZ": "HAMBURG",
    "JA": "JIANGXI",
    "JB": "BRANDENBURG",
    "JC": "BURSA",
    "JD": "EDIRNE",
    "JE": "ANTALYA",
    "JF": "KHOVD",
    "JG": "NEWYORK",
    "JH": "JHB",
    "JI": "JILIN",
    "JJ": "JEJU",
    "JK": "NOVGOROD",
    "JL": "TURKEY",
    "JM": "SUDAN",
    "JN": "AFGHANISTAN",
    "JO": "JOHANNESBURG",
    "JP": "JAPAN",
    "JQ": "NIGERIA",
    "JR": "KENTUCKY",
    "JS": "JIANGSU",
    "JT": "TALCAHUANO",
    "JU": "UTAH",
    "JV": "MASSACHUSETTS",
    "JW": "KANSAS",
    "JX": "JIANGXIDONGHU",
    "JY": "MONTANA",
    "JZ": "IOWA",
    "KA": "KASAULI",
    "KB": "KAGAWA",
    "KC": "KOCHI",
    "KD": "KANAGAWA",
    "KE": "KUMAMOTO-C",
    "KF": "KUM",
    "KG": "KAGOSHIMA",
    "KH": "KHABAROVSK",
    "KI": "KYONGGI",
    "KJ": "KWANGJU",
    "KK": "KYONGBUK",
    "KL": "KHON-KAEN",
    "KM": "KUMAMOTO",
    "KN": "KANGWON",
    "KO": "KOBE",
    "KP": "NAKHON-PATHOM",
    "KQ": "KAMPHAENG-PHET",
    "KR": "KOREA",
    "KS": "KAWASAKI",
    "KT": "KYOTO-C",
    "KU": "KITAKYUSYU",
    "KV": "KANCHANABURI",
    "KW": "KWANGJIN",
    "KX": "KYUNGNAM",
    "KY": "KITAKYUSHU",
    "KZ": "KYONGNAM",
    "LA": "LAUSANNE",
    "LB": "LOPBURI",
    "LC": "LOP-BURI",
    "LD": "LYON-TRS",
    "LE": "LENINGRAD",
    "LF": "LOBBURI",
    "LG": "LIAONING",
    "LH": "LYON-CHU",
    "LI": "LINNKOPING",
    "LJ": "LAMPANG",
    "LK": "LAMOANG",
    "LL": "CASTILLA",
    "LM": "LIMOGES",
    "LN": "LEON",
    "LO": "LOSANGELES",
    "LP": "LIPETZK",
    "LQ": "LOEI",
    "LR": "LA-REUNION",
    "LS": "LISBON",
    "LT": "LATVIA",
    "LU": "LOUISIANA",
    "LV": "SLOVENIA",
    "LW": "SCHLESWIG-HOLSTEIN",
    "LX": "LEICESTERSHIRE",
    "LY": "LYON",
    "LZ": "LINCOLN",
    "MA": "MADRID",
    "MB": "MOROCCO",
    "MC": "MICHIGAN",
    "MD": "MAD",
    "ME": "MEMPHIS",
    "MF": "CLERMONTFERRAND",
    "MG": "MIYAGI",
    "MH": "MAE-HONG-SORN",
    "MI": "MISSISSIPPI",
    "MJ": "MIE",
    "MK": "MIYAZAKI",
    "ML": "MALY",
    "MM": "MALMO",
    "MN": "MAE-HONG-SON",
    "MO": "MISSOURI",
    "MP": "MONTPELLIER",
    "MQ": "MACAU",
    "MR": "MADAGASCAR",
    "MS": "MINNESOTA",
    "MT": "SAMUT-SAKHON",
    "MU": "SAMUT-PRAKAN",
    "MV": "MONGOLIA",
    "MW": "MOSCOW",
    "MX": "MEXICO",
    "MY": "MAYOCLINIC",
    "MZ": "MASSACHUSETS",
    "NA": "NANCHANG",
    "NB": "NEBRASKA",
    "NC": "NEWCASTLE",
    "ND": "NORTH-DAKOTA",
    "NE": "NICE",
    "NF": "NRW",
    "NG": "NIIGATA",
    "NH": "NARA",
    "NI": "NIJMEGEN",
    "NJ": "NEW-JERSEY",
    "NK": "NAGANO",
    "NL": "NETHERLANDS",
    "NM": "NAGOYA",
    "NN": "NINGBO",
    "NO": "NORTH-CAROLINA",
    "NP": "NIGATA",
    "NQ": "NAPAL",
    "NR": "NOVOSIBIRSK",
    "NS": "NAGASAKI",
    "NT": "NIIGATA-C",
    "NU": "NONTHABURI",
    "NV": "NEVADA",
    "NW": "NEW-CALEDONIA",
    "NX": "NINGXIA",
    "NY": "NEW-YORK",
    "NZ": "NAKHON-RATCHASIMA",
    "OA": "OSAKA-C",
    "OB": "NOVY-BYDZOV",
    "OC": "OTAGA",
    "OD": "NORDRHEIN",
    "OE": "ORENSE",
    "OF": "OREL",
    "OG": "OREGON",
    "OH": "OHIO",
    "OI": "OITA",
    "OJ": "BOLIVIA",
    "OK": "OKLAHOMA",
    "OL": "CHOIBALSAN",
    "OM": "OMSK",
    "ON": "OKINAWA",
    "OO": "NAKONNAYOK",
    "OP": "GUADELOUPE",
    "OQ": "CHONGQING",
    "OR": "ORADEA",
    "OS": "OSLO",
    "OT": "OTAGO",
    "OU": "OUJDA",
    "OV": "OVIEDO",
    "OW": "NORWAY",
    "OX": "MAINE",
    "OY": "OKAYAMA",
    "OZ": "RUSSIA",
    "PA": "PARIS",
    "PB": "PARMA",
    "PC": "PORT-CHALMERS",
    "PD": "POLAND",
    "PE": "PERTH",
    "PF": "PHATTHALUNG",
    "PG": "PERUGIA",
    "PH": "PHILIPPINES",
    "PI": "POITIERS",
    "PJ": "PRAJIANBURI",
    "PK": "PHUKET",
    "PL": "PILSEN",
    "PM": "PANAMA",
    "PN": "PUSAN",
    "PO": "PUERTO-RICO",
    "PP": "PRACHUAP-KHIRI-KHAN",
    "PQ": "PATHUMTHANI",
    "PR": "PRAGUE",
    "PS": "PENNSYLVANIA",
    "PT": "PATHUM-THANI",
    "PU": "PERU",
    "PV": "PHILLIPINES",
    "PW": "PRACHINBURI",
    "PX": "PHETCHABUN",
    "PY": "PARAGUAY",
    "PZ": "PRABCHINBURI",
    "QA": "QUANZHOU",
    "QB": "BANGLADESH",
    "QC": "CORDOBA",
    "QD": "ARKANSAS",
    "QE": "MARTINIQUE",
    "QF": "COLOMBIA",
    "QG": "JAMAICA",
    "QH": "POL",
    "QI": "KIEV",
    "QJ": "ODESSA",
    "QK": "PALENCIA",
    "QL": "NOUACKCHOTT",
    "QM": "HESSEN",
    "QN": "GERONA",
    "QO": "MAURITIUS",
    "QP": "FUKUOKAC",
    "QQ": "CONGQING",
    "QR": "GIFUC",
    "QS": "HIROSHIMAC",
    "QT": "QINGDAO",
    "QU": "QUEENSLAND",
    "QV": "KUMAMOTOC",
    "QW": "NIIGATAC",
    "QX": "SENDAIH",
    "QY": "MARSEILLE",
    "QZ": "BORDEAUX",
    "RA": "ROME",
    "RB": "ROI-ET",
    "RC": "PRACHINMURI",
    "RD": "ROTTERDAM",
    "RE": "REUNION",
    "RF": "RABAT",
    "RG": "ARAGON",
    "RH": "RHEINLAND-PFALZ",
    "RI": "SORIA",
    "RJ": "RIO-DE-JANEIRO",
    "RK": "SAMUTPRAKHAN",
    "RL": "BRATISLAVA",
    "RM": "ROMANIA",
    "RN": "SURIN",
    "RO": "ROMA",
    "RP": "PRACHUAPKHIRIKHAN",
    "RQ": "ROSTOV-ON-DON",
    "RR": "SARABURI",
    "RS": "BURGOS",
    "RT": "RATCHABURI",
    "RU": "RU",
    "RV": "ARVAIKHEER",
    "RW": "NARATHIWAT",
    "RX": "ROSTOVDON",
    "RY": "KRASNOYARSK",
    "RZ": "TRABZON",
    "SA": "SOUTH-AUSTRALIA",
    "SB": "SAGA",
    "SC": "SOUTH-CAROLINA",
    "SD": "SHANGDONG",
    "SE": "SENDAI",
    "SF": "SOFIA",
    "SG": "SHIGA",
    "SH": "SHANGHAI",
    "SI": "SICHUAN",
    "SJ": "SPAIN",
    "SK": "ST.-ETIENNE",
    "SL": "SCOTLAND",
    "SM": "SAMARA",
    "SN": "ST-ETIENNE",
    "SO": "SOUTH-DAKOTA",
    "SP": "SINGAPORE",
    "SQ": "SACHSEN-ANHALT",
    "SR": "SAPPORO",
    "SS": "ST.-PETERSBURG",
    "ST": "STOCKHOLM",
    "SU": "SEOUL",
    "SV": "SHANTOU",
    "SW": "SW",
    "SX": "SOPHIA",
    "SY": "SYDNEY",
    "SZ": "SANTIAGO",
    "TA": "TAIWAN",
    "TB": "THURINGEN",
    "TC": "TOCHIGI",
    "TD": "TRENTO",
    "TE": "TEXAS",
    "TF": "TARRAGONA",
    "TG": "TONGA",
    "TH": "THESSALONIKI",
    "TI": "TILBURG",
    "TJ": "TIANJIN",
    "TK": "TAK",
    "TL": "THAILAND",
    "TM": "TASMANIA",
    "TN": "TOULON",
    "TO": "TOULOUSE",
    "TP": "TOYAMA",
    "TQ": "TOKUSHIMA",
    "TR": "TEHRAN",
    "TS": "TRIESTE",
    "TT": "TOTTORI",
    "TU": "TULA",
    "TV": "TOWNSVILLE",
    "TW": "TENNESSEE",
    "TX": "TX",
    "TY": "TOKYO",
    "TZ": "TANGER",
    "UA": "U.K.",
    "UB": "UTTARADIT",
    "UC": "UD",
    "UD": "UDORN",
    "UE": "UBON-RATCHATHANI",
    "UF": "UBONRATCHATHANI",
    "UG": "URUAGUAY",
    "UH": "SUPHANBURI",
    "UI": "UNITED-KINGDOM",
    "UJ": "UTHAI-THANI",
    "UK": "UK",
    "UL": "ULSAN",
    "UM": "UMEA",
    "UN": "UNITEDKINGDOM",
    "UO": "SUKHOTHAI",
    "UP": "SAMUTPRAKAN",
    "UQ": "UMEA",
    "UR": "URUGUAY",
    "US": "USSR/RUSSIA",
    "UT": "UTRECHT",
    "UU": "ULAN-UDE",
    "UV": "SUCEAVA",
    "UW": "ULAANBAATAR",
    "UX": "BENELUX",
    "UY": "GUYANE",
    "UZ": "UKRAINE",
    "VA": "VALENCIA",
    "VB": "VINA-DEL-MAR",
    "VC": "COSTARICA",
    "VD": "VALLADOLID",
    "VE": "VIETNAM",
    "VF": "COOK-ISLAND",
    "VG": "VIRGINIA",
    "VH": "JEONBUK",
    "VI": "VICTORIA",
    "VJ": "GUATEMALA",
    "VK": "SLOVAKIA",
    "VL": "VOLGOGRAD",
    "VM": "VLADIMIR",
    "VN": "VIENNA",
    "VO": "VORONEZH",
    "VP": "KENYA",
    "VQ": "SING",
    "VR": "STAVROPOL",
    "VT": "VERMONT",
    "VV": "AVILA",
    "VY": "IVORY-COAST",
    "VZ": "VENEZUELA",
    "WA": "WASHINGTON",
    "WB": "WY--",
    "WE": "WELLINGTON",
    "WG": "GWANGJU",
    "WH": "WUZHOU",
    "WK": "WAIKATO",
    "WM": "WAKAYAMA",
    "WN": "WISCONSIN",
    "WO": "GANGWON",
    "WR": "WARSAW",
    "WS": "WEST-VIRGINIA",
    "WU": "WUHAN",
    "WV": "WESTVIRGINIA",
    "WW": "NEWMARKET",
    "WX": "NEW-MEXICO",
    "WY": "WYOMING",
    "WZ": "SWITZERLAND",
    "XI": "XINJIANGCHANGJI",
    "XJ": "XINJIANG-HUTUBI",
    "XM": "XIAMEN",
    "XN": "XINJIANGHUTUBI",
    "XX": "CX-ROUSSE",
    "YA": "YAMAGA",
    "YB": "YUNAN",
    "YD": "BYDGOSZCZ",
    "YE": "GYEONGNAM",
    "YG": "YAMAGUCHI",
    "YI": "YAMANESHI",
    "YK": "YOKOSUKA",
    "YL": "MARYLAND",
    "YM": "YAMANASHI",
    "YN": "KYUNGBUK",
    "YO": "YOKOHAMA",
    "YR": "YARYSLAVL",
    "YS": "YAROSLAVL",
    "YT": "YAMAGATA",
    "YU": "YUNNAN",
    "YY": "NAKHONNAYOK",
    "YZ": "YUNGNAM",
    "ZA": "ZAMBIA",
    "ZB": "AZERBAIJAN",
    "ZE": "ZAGREB",
    "ZG": "ZARAGOSSA",
    "ZH": "ZHEJIANG",
    "ZI": "IZMIR",
    "ZL": "ZLIN",
    "ZM": "ZAMORA",
    "ZN": "RYAZAN",
    "ZR": "ZARAGOZA",
    "ZU": "SHIZUOKAC",
    "ZX": "CHUNGBUK",
    "ZY": "SAITAMA",
    "GK": "GOETEBORG"
}
