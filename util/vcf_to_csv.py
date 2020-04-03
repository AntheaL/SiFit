import numpy as np
import vcfpy
import argparse
import time

def get_generators(calls, dl='/'):
    return sorted(list(set(
        [
            int(n) for call in calls if call.data['GT']
            for n in call.data['GT'].split(dl) if n not in {'.', '0'}
        ]
    )))

def pl_to_likelihood(pl, out_sum=True):
    pl = np.array(pl).astype(float)
    x_m = 10*np.log(
        np.sum(np.power(10, -pl/10))
    )/np.log(10)
    pl += x_m
    likelihood = np.power(10, -pl/10)
    if out_sum:
        return np.sum(likelihood[1:])
    return likelihood


def get_mutations(calls, generators, dl='/', use_pl = None):
    if use_pl:
        mutations = np.zeros((1, len(calls))).astype(float)
    else:
        mutations = np.zeros((len(generators), len(calls))).astype(int)
    first = True
    for i, call in enumerate(calls):
        if not call.data.get('AD'):
            mutations[:, i] = 3
            continue
        if use_pl:
            pl = call.data['PL']
            mutations[0, i] = pl_to_likelihood(pl)
            if first:
                print(mutations[0,i])
                first=False

        else:
            gb = call.data['GT'].split(dl)
            for x in gb:
                x = int(x)
                if x:
                    mutations[generators.index(x), i] += 1
    return mutations

def main():

    parser = argparse.ArgumentParser(
        description="Script to retreive mutation matrix for SCITE from vcf file"
    )
    parser.add_argument("--vcf-file", help="input vcf file", required=True, type=str)
    parser.add_argument(
        "--dst-path",
        help="full path where to save output",
        required=False,
        type=str,
        default='mutations.csv'
    )
    parser.add_argument(
        "--filter-thresh",
        help="filter mutations appearing less times than given threshold",
        required=False,
        type=int,
    )
    parser.add_argument(
        "--use-pl",
        action="store_true",
        required=False
    )
    parser.add_argument(
        "--differentiated",
        help="not to store mutations in a hierarchical way",
        required=False,
        action='store_true'
    )
    args = parser.parse_args()

    if args.differentiated:
        print("Will not store mutations in a hierarchical order.")

    reader = vcfpy.Reader.from_path(args.vcf_file)
    mutations = None
    n_per_site = []

    t0 = time.time()
    i=0
    for record in reader:
        calls = record.calls
        generators = get_generators(calls)
        if i==0:
            for call in calls:
                if call.data.get('AD'):
                    print(call.data)
                    break
        curr_mutations = get_mutations(calls, generators, use_pl = args.use_pl)
        n_per_site.append(curr_mutations.shape[0])
        if mutations is None:
            mutations = curr_mutations
        else:
            mutations = np.concatenate([mutations, curr_mutations], axis=0)
        if i==0:
            print(f"time taken per record: {time.time()-t0}")
        if i==10:
            print(f"current shape: {mutations.shape}")
            print(f"first elements: {mutations[0,:5]}")
        i+=1

    if args.filter_thresh:
        mbis = mutations.copy()
        for x in mbis:
            x[x==3] = 0
        fltr = np.where([np.sum(x)>args.filter_thresh for x in mbis])[0]
        mutations = mutations[fltr]

    print(f"fineal shape of mutation matrix: {mutations.shape}")

    np.savetxt(args.dst_path, mutations, delimiter=' ', fmt='%.4f')
    #np.save("np_" + args.dst_path, mutations)

    print(f"proportion of missing data: {len(mutations[mutations==3])/np.prod(mutations.shape)}")

if __name__ == "__main__":
    main()
