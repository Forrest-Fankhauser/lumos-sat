"""
Tools for getting TLEs from Celestrak and SpaceTrack
"""

import requests

def query_celestrak(filename, query, value, format, supplemental = False):
    """
    Tool to interface with Celestrak API.
    See `celestrak.org <https://celestrak.org/NORAD/documentation/gp-data-formats.php>`_ for more details.

    :param filename: File where downloaded data is written
    :type filename: str
    :param query: CATNR or INTDES or GROUP or NAME or SPECIAL
    :type query: str
    :param value: Determines how the query is filtered
    :type value: str
    :param format: TLE or 2LE or XML or KVN or JSON or JSON-PRETTY or CSV
    :type format: str
    :param supplemental: Whether to pull from Celestrak's suplemental GP data
    :type supplemental: bool, optional

    """

    url = (f"https://celestrak.org/NORAD/elements/"
           f"{'supplemental/sup-' if supplemental else ''}"
           f"gp.php?{query}={value}&FORMAT={format}")
    
    response = requests.get(url)

    with open(filename, "wb") as file:
        file.write(response.content)