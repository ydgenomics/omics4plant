#!/usr/bin/env python3
# from: yinzhanhao
import gzip, os, sys, multiprocessing as mp

def stream(in_path):
    with gzip.open(in_path, 'rt') as f:
        while True:
            h = f.readline().rstrip('\n\r')
            if not h: break
            s = f.readline().rstrip('\n\r')
            p = f.readline().rstrip('\n\r')
            q = f.readline().rstrip('\n\r')
            if not q: break
            yield h, s, p, q

# 真正裁剪函数
def _trim(args):
    in_path, out_path, seqlen, regions = args
    tot = skip = 0
    with gzip.open(out_path, 'wt') as out:
        for h, s, p, q in stream(in_path):
            tot += 1
            if len(s) < seqlen:
                skip += 1
                continue
            new_s = ''.join(s[start:end] for start, end in regions)
            new_q = ''.join(q[start:end] for start, end in regions)
            out.write(f'{h}\n{new_s}\n{p}\n{new_q}\n')
    return tot, skip

def main(r1_in, r2_in, out_dir):
    r1_out = os.path.join(out_dir, os.path.basename(r1_in).replace('_1.fq.gz', '_trim_1.fq.gz'))
    r2_out = os.path.join(out_dir, os.path.basename(r2_in).replace('_2.fq.gz', '_trim_2.fq.gz'))

    # 参数包：(输入, 输出, 最小长度, 截取区间)
    r1_job = (r1_in, r1_out, 20, [(0, 20)])
    r2_job = (r2_in, r2_out, 42, [(0, 10), (16, 26), (32, 42)])

    with mp.Pool(2) as pool:
        (tot1, skip1), (tot2, skip2) = pool.map(_trim, [r1_job, r2_job])

    print(f'[INFO] R1 total={tot1} skip={skip1}  R2 total={tot2} skip={skip2}  ok={min(tot1-skip1, tot2-skip2)}')

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])