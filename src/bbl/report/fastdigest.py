"""Thermo Scientific FastDigest reaction conditions, one row per enzyme.

Scraped from the supplier's own table on 2026-08-15:
https://www.thermofisher.com/us/en/home/brands/thermo-scientific/molecular-biology/thermo-scientific-restriction-modifying-enzymes/restriction-enzymes-thermo-scientific/fastdigest-thermo-scientific/reaction-conditions-for-fastdigest-enzymes.html

This is a **cache of a supplier document, not knowledge**. It goes stale when Thermo
reformulates an enzyme, so every table rendered from it carries the retrieval date and the
link. An enzyme that is not in here renders as blank cells rather than a guess -- there are
176 FastDigest enzymes and rather more restriction enzymes in the world.

Row: ``(reaction temperature C, minutes for 1 ug plasmid in 20 uL, bp from the end of a
fragment needed for complete digestion, thermal inactivation, hours before star activity)``.
"""

from __future__ import annotations

#: Date the table above was retrieved. Rendered next to the numbers.
RETRIEVED = "2026-08-15"

SOURCE_URL = (
    "https://www.thermofisher.com/us/en/home/brands/thermo-scientific/molecular-biology/"
    "thermo-scientific-restriction-modifying-enzymes/restriction-enzymes-thermo-scientific/"
    "fastdigest-thermo-scientific/reaction-conditions-for-fastdigest-enzymes.html"
)

#: enzyme -> (temp_c, plasmid_minutes, bp_from_end, thermal_inactivation, star_free_hours)
CONDITIONS: dict[str, tuple[str, str, str, str, str]] = {
    'AanI': ('37', '5', '5', '65°C, 5 min', '16'),
    'AasI': ('37', '5', '1', '80°C, 10 min', '1'),
    'AatII': ('37', '20', '5', '80°C, 5 min', '16'),
    'Acc65I': ('37', '5', '3', '65°C, 5 min', '16'),
    'AdeI': ('37', '5', '2', '80°C, 5 min', '1'),
    'AjuI': ('37', '5', 'nd', '65°C, 5 min', '16'),
    'AluI': ('37', '15', '4', '65°C, 5 min', '16'),
    'Alw21I': ('37', '5', '1', '80°C, 20 min', '16'),
    'Alw26I': ('37', '5', '1', '65°C, 5 min', '16'),
    'Alw44I': ('37', '60', '4', '80°C, 5 min', '16'),
    'ApaI': ('37', '5', '2', '65°C, 5 min', '16'),
    'BamHI': ('37', '5', '2', '80°C, 5 min', '1'),
    'BclI': ('37', '5', '3', '80°C, 20 min', '16'),
    'BcnI': ('37', '5', '3', '80°C, 20 min', '16'),
    'BcuI': ('37', '5', '1', 'No (chloroform extraction)', '16'),
    'BfmI': ('37', '5', '2', '65°C, 20 min', '16'),
    'BfoI': ('37', '5', '2', '65°C, 10 min', '16'),
    'BglI': ('37', '5', '3', '65°C, 5 min', '2'),
    'BglII': ('37', '20', '3', 'No (chloroform extraction)', '16'),
    'Bme1390I': ('37', '5', '4', 'No (chloroform extraction)', '16'),
    'BmsI': ('37', '5', '1', '65°C, 5 min', '16'),
    'BoxI': ('37', '5', '3', '80°C, 20 min', '16'),
    'BpiI': ('37', '5', '1', '65°C, 10 min', '16'),
    'BplI': ('37', '20', 'nd', '65°C, 5 min', '16'),
    'Bpu10I': ('37', '15', '3', '80°C, 5 min', '1'),
    'Bpu1102I': ('37', '5', '2', '80°C, 5 min', '16'),
    'BseDI': ('37', '5', '3', '80°C, 5 min', '16'),
    'BseGI': ('37', '5', '2', '80°C, 5 min', '16'),
    'BseJI': ('65', '5', '2', 'No (chloroform extraction)', '1'),
    'BseLI': ('37', '5', '2', 'No (chloroform extraction)', '16'),
    'BseMI': ('37', '15', '3', '80°C, 5 min', '16'),
    'BseMII': ('55', '30', '2', '80°C, 5 min', '16'),
    'BseNI': ('65', '5', '2', '80°C, 5 min', '16'),
    'BseSI': ('37', '5', '3', '80°C, 15 min', '16'),
    'BseXI': ('65', '15', '3', '80°C, 15 min', '16'),
    'Bsh1236I': ('37', '5', '1', '80°C, 10 min', '16'),
    'Bsh1285I': ('37', '15', '3', '80°C, 15 min', '1'),
    'BshNI': ('37', '60', '2', '65°C, 10 min', '16'),
    'BshTI': ('37', '5', '3', '80°C, 5 min', '16'),
    'Bsp119I': ('37', '5', '1', '80°C, 5 min', '16'),
    'Bsp120I': ('37', '5', '3', '80°C, 10 min', '16'),
    'Bsp1407I': ('37', '5', '3', '80°C, 10 min', '16'),
    'Bsp143I': ('37', '5', '4', '65°C, 20 min', '16'),
    'BspLI': ('37', '5', '2', '65°C, 20 min', '16'),
    'BspOI': ('37', '5', '3', '80°C, 5 min', '1'),
    'BspTI': ('37', '5', '4', 'No (chloroform extraction)', '16'),
    'Bst1107I': ('37', '5', '3', 'No (chloroform extraction)', '6'),
    'BstXI': ('37', '5', '4', '80°C, 5 min', '1'),
    'Bsu15I': ('37', '5', '2', '65°C, 15 min', '16'),
    'BsuRI': ('37', '5', '3', 'No (chloroform extraction)', '16'),
    'BveI': ('37', '15', '5', '65°C, 5 min', '16'),
    'CaiI': ('37', '5', '4', '65°C, 10 min', '16'),
    'Cfr10I': ('37', '5', '2', 'No (chloroform extraction)', '2'),
    'Cfr13I': ('37', '5', '2', 'No (chloroform extraction)', '16'),
    'CpoI': ('37', '5', '1', '65°C, 5 min', '16'),
    'CseI': ('37', '15', '2', '80°C, 5 min', '16'),
    'CsiI': ('37', '10', '4', '65°C, 5 min', '16'),
    'Csp6I': ('37', '5', '2', '80°C, 10 min', '16'),
    'DpnI': ('37', '5', '1', '80°C, 5 min', '16'),
    'DraI': ('37', '5', '4', '65°C, 5 min', '16'),
    'Eam1104I': ('37', '5', '3', '80°C, 5 min', '6'),
    'Eam1105I': ('37', '5', '3', '65°C, 5 min', '16'),
    'Ecl136II': ('37', '5', '1', '65°C, 5 min', '6'),
    'Eco105I': ('37', '5', '3', '65°C, 5 min', '16'),
    'Eco130I': ('37', '5', '3', '65°C, 5 min', '16'),
    'Eco147I': ('37', '5', '3', '80°C, 10 min', '16'),
    'Eco31I': ('37', '5', '3', '65°C, 5 min', '16'),
    'Eco32I': ('37', '5', '2', 'No (chloroform extraction)', '16'),
    'Eco47I': ('37', '5', '1', '80°C, 20 min', '16'),
    'Eco47III': ('37', '5', '4', '65°C, 5 min', '16'),
    'Eco52I': ('37', '20', '3', '65°C, 5 min', '16'),
    'Eco57I': ('37', '30', '2', '65°C, 5 min', '16'),
    'Eco72I': ('37', '5', '2', '80°C, 10 min', '1'),
    'Eco81I': ('37', '5', '2', '80°C, 10 min', '16'),
    'Eco88I': ('37', '5', '2', '65°C, 5 min', '16'),
    'Eco91I': ('37', '5', '2', '65°C, 10 min', '2'),
    'EcoO109I': ('37', '10', '2', '65°C, 5 min', '16'),
    'EcoRI': ('37', '5', '2', '80°C, 5min', '0.5'),
    'EheI': ('37', '5', '2', '65°C, 5 min', '6'),
    'Esp3I': ('37', '5', '1', '65°C, 10 min', '6'),
    'FaqI': ('37', '30', '3', '65°C, 5 min', '16'),
    'FokI': ('37', '5', '1', '65°C, 5 min', '1'),
    'FspAI': ('37', '5', '4', '65°C, 5 min', '16'),
    'FspBI': ('37', '15', '3', '80°C, 5 min', '16'),
    'GsuI': ('30', '30', '2', '65°C, 5 min', '16'),
    'HhaI': ('37', '5', '1', 'No (chloroform extraction)', '16'),
    'Hin1I': ('37', '5', '1', '65°C, 5 min', '16'),
    'Hin1II': ('37', '10', '4', '80°C, 5 min', '16'),
    'Hin6I': ('37', '10', '2', '80°C, 10 min', '16'),
    'HincII': ('37', '5', '1', '65°C, 5 min', '16'),
    'HindIII': ('37', '5', '3', '80°C, 10 min', '16'),
    'HinfI': ('37', '5', '2', '65°C, 20 min', '16'),
    'HpaII': ('37', '5', '3', '65°C, 5 min', '16'),
    'Hpy8I': ('37', '5', '2', '80°C, 5 min', '16'),
    'HpyF10VI': ('37', '5', '2', '80°C, 5 min', '16'),
    'HpyF3I': ('37', '5', '3', '65°C, 5 min', '16'),
    'KflI': ('37', '5', '5', 'No (chloroform extraction)', '16'),
    'Kpn2I': ('37', '5', '2', '80°C, 5 min', '16'),
    'KpnI': ('37', '5', '3', '80°C, 5 min', '16'),
    'KspAI': ('37', '5', '3', '65°C, 20 min', '0.5'),
    'LguI': ('37', '5', '2', '65°C, 5 min', '16'),
    'Lsp1109I': ('37', '5', '2', '65°C, 5 min', '1'),
    'MauBI': ('37', '5', '5', '65°C, 5 min', '16'),
    'MbiI': ('37', '5', '3', '65°C, 5 min', '16'),
    'MboI': ('37', '5', '1', '65°C, 15 min', '16'),
    'MboII': ('37', '5', '2', '65°C, 5 min', '16'),
    'MlsI': ('37', '5', '3', '80°C, 20 min', '16'),
    'MluI': ('37', '5', '3', '80°C, 5 min', '16'),
    'MnlI': ('37', '5', '2', '65°C, 5 min', '16'),
    'Mph1103I': ('37', '5', '3', '65°C, 15 min', '6'),
    'MreI': ('37', '5', '3', '80°C, 5 min', '16'),
    'MspI': ('37', '5', '3', 'No (chloroform extraction)', '16'),
    'MssI': ('37', '5', '2', '65°C, 10 min', '16'),
    'MunI': ('37', '5', '3', 'No (chloroform extraction)', '16'),
    'Mva1269I': ('37', '5', '3', '65°C, 5 min', '16'),
    'MvaI': ('37', '5', '4', 'No (chloroform extraction)', '1'),
    'NcoI': ('37', '10', '3', '65°C, 15 min', '16'),
    'NdeI': ('37', '5', '3', '65°C, 5 min', '6'),
    'NheI': ('37', '15', '5', '65°C, 5 min', '6'),
    'NmuCI': ('37', '5', '2', '65°C, 5 min', '16'),
    'NotI': ('37', '30', '2', '80°C, 5 min', '16'),
    'NsbI': ('37', '5', '3', '65°C, 15 min', '16'),
    'OliI': ('37', '5', '3', '65°C, 5 min', '6'),
    'PacI': ('37', '5', '2', '65°C, 10 min', '16'),
    'PaeI': ('37', '5', '5', '65°C, 5 min', '16'),
    'PagI': ('37', '10', '3', '80°C, 5 min', '0.5'),
    'PdiI': ('37', '5', '4', '65°C, 15 min', '16'),
    'PdmI': ('37', '5', '2', '65°C, 5 min', '16'),
    'PfeI': ('37', '5', '2', '65°C, 5 min', '16'),
    'Pfl23II': ('37', '15', '3', '65°C, 5 min', '16'),
    'PfoI': ('37', '5', '3', '65°C, 5 min', '16'),
    'Ppu21I': ('37', '5', '1', '65°C, 5 min', '16'),
    'Psp1406I': ('37', '5', '3', '65°C, 5 min', '16'),
    'Psp5II': ('37', '5', '2', '80°C, 5 min', '6'),
    'PspFI': ('37', '5', '4', '80°C, 5 min', '16'),
    'PstI': ('37', '5', '3', 'No (chloroform extraction)', '16'),
    'PsuI': ('37', '5', '2', '80°C, 5 min', '16'),
    'PsyI': ('37', '5', '2', '80°C, 5 min', '4'),
    'PteI': ('37', '5', '3', '80°C, 5 min', '16'),
    'PvuI': ('37', '15', '4', '80°C, 5 min', '16'),
    'PvuII': ('37', '5', '1', 'No (chloroform extraction)', '0.5'),
    'RruI': ('37', '5', '5', 'No (chloroform extraction)', '4'),
    'RsaI': ('37', '5', '3', 'No (chloroform extraction)', '16'),
    'RseI': ('37', '5', '4', '65°C, 20 min', '16'),
    'SacI': ('37', '15', '1', '65°C, 5 min', '16'),
    'SalI': ('37', '5', '3', '65°C, 10 min', '16'),
    'SaqAI': ('37', '5', '4', '65°C, 5 min', '1'),
    'SatI': ('37', '5', '3', '65°C, 5 min', '16'),
    'ScaI': ('37', '5', '4', '65°C, 10 min', '16'),
    'SchI': ('37', '5', '2', '80°C, 5 min', '1'),
    'SdaI': ('37', '5', '3', 'No (chloroform extraction)', '1'),
    'SduI': ('37', '5', '2', '80°C, 15 min', '1'),
    'SfaAI': ('37', '5', '5', '80°C, 5 min', '16'),
    'SfiI': ('50', '15', '5', 'No (chloroform extraction)', '16'),
    'SgsI': ('37', '5', '2', '65°C, 20 min', '16'),
    'SmaI': ('37', '5', '1', '65°C, 5 min', '16'),
    'SmiI': ('37', '20', '2', '65°C, 15 min', '1'),
    'SsiI': ('37', '5', '2', '65°C, 5 min', '4'),
    'SspI': ('37', '5', '2', '65°C, 5 min', '1'),
    'TaaI': ('65', '5', '3', 'No (chloroform extraction)', '16'),
    'TaiI': ('65', '5', '2', 'No (chloroform extraction)', '16'),
    'TaqI': ('65', '5', '4', 'No (chloroform extraction)', '16'),
    'TasI': ('65', '5', '3', 'No (chloroform extraction)', '6'),
    'TatI': ('65', '15', '5', 'No (chloroform extraction)', '1'),
    'TauI': ('55', '60', '3', 'No (chloroform extraction)', '16'),
    'Tru1I': ('65', '5', '3', 'No (chloroform extraction)', '2'),
    'TscAI': ('65', '5', '2', 'No (chloroform extraction)', '6'),
    'Van91I': ('37', '5', '3', '65°C, 10 min', '6'),
    'VspI': ('37', '5', '2', '65°C, 5 min', '16'),
    'XagI': ('37', '5', '2', '65°C, 5 min', '16'),
    'XapI': ('37', '5', '4', '80°C, 5 min', '16'),
    'XbaI': ('37', '5', '2', '65°C, 20 min', '16'),
    'XceI': ('37', '5', '2', '65°C, 5 min', '16'),
    'XhoI': ('37', '5', '2', '80°C, 5 min', '16'),
    'XmaJI': ('37', '10', '2', 'No (chloroform extraction)', '16'),
    'XmiI': ('37', '5', '2', '65°C, 5 min', '16'),
}


def lookup(enzyme: str):
    """Conditions for ``enzyme``, trying its isoschizomers before giving up.

    Our planner names enzymes by the NEB-preferred label (MfeI, ClaI); FastDigest sells the
    same activity under another name (MunI, Bsu15I). Biopython knows the synonym set, so the
    mapping is derived rather than maintained by hand.

    Returns ``(fastdigest_name, row)`` or ``(None, None)``.
    """
    name = str(enzyme)
    if name in CONDITIONS:
        return name, CONDITIONS[name]
    try:
        from Bio import Restriction

        for other in getattr(Restriction, name).isoschizomers():
            if str(other) in CONDITIONS:
                return str(other), CONDITIONS[str(other)]
    except (AttributeError, ImportError):
        pass
    return None, None
