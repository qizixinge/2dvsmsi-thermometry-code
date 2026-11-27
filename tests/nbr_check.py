"""Independent ROI neighbourhood recomputation from voxel maps."""
import numpy as np


def percentile(values, center):
    less = np.count_nonzero(values < center)
    ties = np.count_nonzero(values == center)
    return 100 * (less + .5 * (ties - 1)) / (values.size - 1)


def agree(actual, expected, hint):
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-9,
                               equal_nan=False, err_msg=hint)


def common_rms(series, temperatures, indices):
    # series[roi, window, time]. Each fit withholds the same time for all ROIs.
    total = np.zeros(indices[0].shape)
    for t in range(9):
        train = np.arange(9) != t
        x = temperatures[train]
        xc = x - x.mean(axis=0)
        y = series[:, :, train]
        yc = y - y.mean(axis=2, keepdims=True)
        cross = (yc * xc.T[:, None, :]).sum(axis=2)
        slope = sum(cross[k, indices[k]] for k in range(4)) / np.sum(xc**2)
        delta = series[:, :, t] - y.mean(axis=2)
        for k in range(4):
            error = delta[k, indices[k]] / slope + x[:, k].mean() - temperatures[t, k]
            total += error**2
    return np.sqrt(total / 36)


def nbr_values(s, r, sos, mip, prf_signal, loo_roi):
    # From verify_statistics.py pass -dphase/(B0*2*np.pi*GAMMA_BAR*1e6*dTE).
    names = ('Gre', 'PfaSos', 'FullSos', 'PfaMip', 'FullMip')
    temps = [r[k] for k in ('temprp1', 'temp1', 'temp3', 'temp1', 'temp3')]
    centers = [s[k] for k in ('CPRF', 'C2D', 'C2D_full', 'C2D', 'C2D_full')]
    maps = (prf_signal, sos[..., 0], sos[..., 2], mip[..., 0], mip[..., 2])
    archived = np.stack([s['ALPHA_PAPER'] * 1e8] +
                        [s[k] for k in ('m2_p_sos','m2_f_sos','m2_p_mip','m2_f_mip')])
    slopes = np.full((5, 4, 5, 5), np.nan)
    rmse = np.full_like(slopes, np.nan)
    valid = np.zeros(slopes.shape, bool)
    signals = np.full(slopes.shape + (9,), np.nan)
    for v, vol in enumerate(maps):
        top, bottom, left, right = s['geometry']['rod_128' if v == 0 else 'rod_100']
        for k, (row, col) in enumerate(np.asarray(centers[v], int)):
            x = temps[v][:, k]
            dx = x - x.mean()
            for i, j in np.ndindex(5, 5):
                a, b = row+i-2, col+j-2  # one-based shifted centre
                if a < 2 or b < 2 or a >= vol.shape[0] or b >= vol.shape[1]:
                    continue
                if a+1 >= top and a-1 <= bottom and b+1 >= left and b-1 <= right:
                    continue
                block = vol[a-2:a+1, b-2:b+1]
                if v and not np.all(np.isfinite(block) & (block > 0)):
                    continue
                y = block.mean(axis=(0, 1), dtype=float)
                signals[v,k,i,j] = y
                slopes[v,k,i,j] = dx @ (y-y.mean()) / (dx @ dx) * (1e8 if v == 0 else 1)
                errors = loo_roi(y[:, None], x[:, None])
                rmse[v,k,i,j] = np.sqrt(np.mean(errors**2))
                valid[v,k,i,j] = True
    assert valid[:,:,2,2].all(), 'Check archived centres, rod box and map validity'
    agree(slopes[:,:,2,2], archived, 'Check map order, GRE sign/units and ROI centres')
    agree(rmse[:,:,2,2], s['loocv_rmse'], 'Check leave-one-time-out fits and probe assignment')
    inner = valid[:,:,1:4,1:4]
    counts = inner[1:].sum(axis=(2,3))
    assert np.all(counts == 9), 'Check radius-1 bounds, rod box and map-specific mask'
    combinations = counts.prod(axis=1)
    assert np.all(combinations == combinations[0]), 'Check per-variant Cartesian products'
    E = dict(NbrNumWin=str(counts[0,0]), NbrNumComb=str(combinations[0]))
    A = dict(slope=slopes, rmse=rmse, valid=valid, cv=[], rms=[], crms=[])
    assert np.all(inner[0].sum(axis=(1,2)) == 9), 'GRE: not every radius-1 window kept; report the kept counts instead of "of nine"'
    for k, label in enumerate(('I','II','III','IV')):
        vals = rmse[0,k,1:4,1:4][inner[0,k]]
        E['NbrGreRmseMed'+label] = '%.2f' % np.median(vals)
        E['NbrGreRank'+label] = '%d' % (1 + np.count_nonzero(vals < rmse[0,k,2,2]))   # 1 = lowest RMSE

    for radius in (1, 2):
        section = slice(2-radius, 3+radius)
        width = (2*radius+1)**2
        indices = np.ogrid[tuple(slice(width) for _ in range(4))]
        cv_grids = []
        kept_grids = []
        for v in range(1, 5):
            def grid(array):
                return array[v,:,section,section].reshape(4, width)
            sl = grid(slopes)
            values = [sl[k, indices[k]] for k in range(4)]
            mean = sum(values)/4
            cv = 100*np.sqrt(sum((z-mean)**2 for z in values)/3)/np.abs(mean)
            cv_grids.append(cv)
            kept = np.logical_and.reduce(np.broadcast_arrays(
                *[grid(valid)[k, indices[k]] for k in range(4)]))
            kept_grids.append(kept)
            center = (width//2,)*4
            if radius == 2:
                E['NbrCvMedTwo'+names[v]] = '%.1f' % np.median(cv[kept])
                continue
            idx = np.broadcast_arrays(*indices)
            rms = np.sqrt(sum(grid(rmse)[k, idx[k]]**2 for k in range(4))/4)
            series = signals[v,:,section,section].reshape(4,width,9)
            crms = common_rms(series, temps[v], idx)
            agree(cv[center], 100*archived[v].std(ddof=1)/abs(archived[v].mean()),
                  'Check sample SD and archived combination indices')
            agree(crms[center], np.sqrt(np.mean(s['common_rmse'][v-1]**2)),
                  'Check shared slope, training means and simultaneous holdout')
            for key, prefix, metric, fmt in [('cv','Cv',cv,'%.1f'), ('rms','Rms',rms,'%.2f'),
                                             ('crms','CommonRms',crms,'%.2f')]:
                vector = metric[kept]
                A[key].append(vector)
                E['Nbr'+prefix+'Med'+names[v]] = fmt % np.median(vector)
                E['Nbr'+prefix+'Pct'+names[v]] = '%.0f' % percentile(vector, metric[center])
        for k, label in enumerate(('Pfa','Full')):
            a, b = cv_grids[k], cv_grids[k+2]
            both = kept_grids[k] & kept_grids[k+2]        # windows kept in both SOS and MIP
            assert np.all(np.isfinite(a[both]) & np.isfinite(b[both])), 'CV undefined for a kept combination'
            suffix = 'Two' if radius == 2 else ''
            if radius == 2:
                E['NbrNumCombTwo'+label] = str(both.sum())
            E['NbrSosLtMip'+suffix+label] = '%.0f' % (100*np.mean(a[both] < b[both]))
    return E, A


def nbr_compare(A, NBR, close):
    def check(name, a, b):
        a, b = np.asarray(a), np.asarray(b)
        assert a.shape == b.shape, (name, a.shape, b.shape)
        close(name+' NaN mask', np.isnan(a).astype(int), np.isnan(b).astype(int))
        keep = ~np.isnan(a)
        x, y = a[keep].astype(float), b[keep].astype(float)
        # verify's close() allows 1e-9 absolute + 1e-9 relative and records the raw absolute difference.  A few
        # shifted windows have slopes near zero, so their LOOCV RMSE is large (up to ~1.6e3 degC) and differs
        # between MATLAB and numpy by ~4e-8 degC absolute (~2e-11 relative) (2026-11-27).
        close(name, x, y, atol=1e-9)
    assert np.array_equal(np.ravel(NBR['off']), np.arange(-2, 3)), 'NBR.off must be -2:2 (axis labels of the offset dimensions)'
    for key in ('slope','rmse','valid'):
        check('NBR.'+key, A[key], NBR[key])
    for key in ('cv','rms','crms'):
        assert len(NBR[key]) == 4, 'Expected four MATLAB cell entries for '+key
        for k, (a,b) in enumerate(zip(A[key], NBR[key])):
            check(f'NBR.{key}[{k}]', np.sort(a), np.sort(np.asarray(b,float).ravel()))
