#! /usr/bin/env python

import eval_one
import traceback

if __name__ == '__main__':

    proposal_file = raw_input('Enter path to proposal list:\n')
    with open(proposal_file) as f:
        proposal_list = f.readlines()
    proposal_list = ['{}.apt'.format(item.strip()) for item in proposal_list]

    for proposal in proposal_list:
        try:
            eval_one.doit(proposal)
        except Exception as error:
            trace = traceback.format_exc()
            with open('error_log.dat', 'a') as error_file:
                error_file.write('{}{}\n\n'.format(str(trace), proposal))
