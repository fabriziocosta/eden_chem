#!/usr/bin/env python
"""Provides io for molecules."""

import os
import requests
import logging
logger = logging.getLogger(__name__)


root_uri = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/'
_pubchem_dir_ = 'PUBCHEM'


def get_assay_description(assay_id):
    """get_assay_description."""
    fname = 'AID%s_info.txt' % assay_id
    full_fname = os.path.join(_pubchem_dir_, fname)
    if not os.path.isfile(full_fname):
        query = root_uri
        query += 'assay/aid/%s/summary/JSON' % assay_id
        reply = requests.get(query)
        text = reply.json()['AssaySummaries']['AssaySummary'][0]['Name']
        with open(full_fname, 'w') as file_handle:
            file_handle.write(text)
    else:
        with open(full_fname, 'r') as file_handle:
            text = ''
            for line in file_handle:
                text += line
    return text


def _get_compounds(fname, active, aid, stepsize=50):
    with open(fname, 'w') as file_handle:
        index_start = 0

        reply = requests.get(_make_rest_query(aid, active=active))
        listkey = reply.json()['IdentifierList']['ListKey']
        size = reply.json()['IdentifierList']['Size']

        for chunk, index_end in enumerate(range(0, size + stepsize, stepsize)):
            if index_end is not 0:
                repeat = True
                while repeat:
                    t = 'Chunk %s) Processing compounds %s to %s (%s)' % \
                        (chunk, index_start, index_end - 1, size)
                    logger.debug(t)
                    query = root_uri
                    query += 'compound/listkey/' + str(listkey)
                    query += '/SDF?&listkey_start=' + str(index_start)
                    query += '&listkey_count=' + str(stepsize)
                    reply = requests.get(query)
                    if 'PUGREST.Timeout' in reply.text:
                        logger.debug("PUGREST TIMEOUT")
                    elif "PUGREST.BadRequest" in reply.text:
                        # bad request means that the server throw our
                        # molecule list away.
                        # we just make a new one
                        logger.debug('bad request %s %d %d %d' % (
                            query, chunk, index_end, size))
                        reply = requests.get(
                            _make_rest_query(aid, active=active))
                        listkey = reply.json()['IdentifierList']['ListKey']
                    elif reply.status_code != 200:
                        logger.debug("UNKNOWN ERRA " + query)
                        logger.debug(reply.status_code)
                        logger.debug(reply.text)
                        raise Exception('UNKNOWN')
                    else:  # everything is OK
                        repeat = False
                        file_handle.write(reply.text)

            index_start = index_end

        logger.debug('compounds available in file: ', fname)


def _make_rest_query(assay_id, active=True):
    if active:
        mode = 'active'
    else:
        mode = 'inactive'
    core = 'assay/aid/%s/cids/JSON?cids_type=%s&list_return=listkey' % \
        (assay_id, mode)
    rest_query = root_uri + core
    return rest_query


def _query_db(assay_id, fname=None, active=True, stepsize=50):
    _get_compounds(fname=fname + ".tmp",
                   active=active,
                   aid=assay_id,
                   stepsize=stepsize)
    os.rename(fname + '.tmp', fname)


def download(assay_id, active=True, stepsize=50):
    """download."""
    if not os.path.exists(_pubchem_dir_):
        os.mkdir(_pubchem_dir_)
    if active:
        fname = 'AID%s_active.sdf' % assay_id
    else:
        fname = 'AID%s_inactive.sdf' % assay_id
    full_fname = os.path.join(_pubchem_dir_, fname)
    if not os.path.isfile(full_fname):
        logger.debug('Querying PubChem for AID: %s' % assay_id)
        _query_db(assay_id,
                  fname=full_fname,
                  active=active,
                  stepsize=stepsize)
    else:
        logger.debug('Reading from file: %s' % full_fname)
    return full_fname
