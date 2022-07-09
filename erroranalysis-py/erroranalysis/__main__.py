#! /usr/bin/env python3
from . import *
import argparse
import sys
import os

def create_parser():
    parser = argparse.ArgumentParser(description='Plots the average and error of an estimator as a function of Monte Carlo bins (optionally performs binning analysis).',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    parser.add_argument('--skip', type=int, default=0,
                        help='Skip <skip> number of bins before averaging')
    parser.add_argument('--estimator', type=str, default="",
                        help='Name of estimator to average as it appears in the header of <filename> or <custom_header>')
    parser.add_argument('--custom_header', type=str, default="",
                        help='A comma separated header for <filename> if one does not exist in the file with the same number of headings as columns of data')
    parser.add_argument('--skip_header', type=int, default=0,
                        help='Skip <skip> number of lines in <filename> before attempting to read the file')
    parser.add_argument('--dtype', type=str, default=None,
                        help='`dtype` passed to numpy.genfromtxt()')
    parser.add_argument('--deletechars', type=str, default=" !#$%&'()*+, -./:;<=>?@[\]^{|}~",
                        help='`deletechars` passed to numpy.genfromtxt()')
    parser.add_argument('--pretty_estimator_name', type=str, default="",
                        help='Name of <estimator> formatted for matplotlib plot labels')
    parser.add_argument('--legend', action='store_true',
                        help='Enable legend on plots')
    parser.add_argument('--labels', type=str, default="",
                        help='A comma separated list with labels for each file in <filename> (need `--legend` enabled) otherwise will use filenames as labels')
    parser.add_argument('--savefig', type=str, default="",
                        help='Name for saving plots. Will prepend extension with plot type and <average_method>. i.e. `--savefig estimator.png` will save a binning analysis plot using jackknife as `estimator_binning-analysis_jackknife.png`')
    parser.add_argument('--mplstylefile', type=str, default="default",
                        help='Location of stylefile to use with plotting')
    parser.add_argument('--dpi', type=float, default=None,
                        help='DPI for plots.')

    parser.add_argument('--bin_size', type=int, default=1,
                        help='Size of bin to average data before processing (3x autocorrelation time is a good size)')
    parser.add_argument('--binning_analysis', action='store_true',
                        help='Perform binning analysis.')
    parser.add_argument('--average_method', type=str, default="jackknife",
                        help='Method for performing averages. Choose from [jackknife, bootstrap, running, cumulative]')
    parser.add_argument('--running_window', type=int, default=1,
                        help='Size of window to use with running average (`--average_method running`)')
    parser.add_argument('--bootstrap_size', type=int, default=100,
                        help='Number of resampled datasets to use with bootstrap method (`--average_method bootstrap`)')

    parser.add_argument('filename', type=str, nargs='+',
                        help='File[s] to process.')
    return parser

def main(argv=None):

    """
    :desc: Plots the average, fluctuations, skew, and kurtosis with error for an estimator as a function of Monte Carlo bins (optionally performs binning analysis).
    """
    if argv is None:
        argv = sys.argv

    parser = create_parser()
    args = parser.parse_args(argv[1:])

    if args.labels:
        labels = [l.strip() for l in args.labels.split(",")]
        if len(labels) != len(args.filename):
            msg = "Length of labels does not match length of filenames"
            raise ValueError(msg)

    if args.custom_header:
        custom_header = [l.strip() for l in args.custom_header.split(",")]
    else:
        custom_header = True

    if args.pretty_estimator_name:
        pname = args.pretty_estimator_name
    else:
        pname = args.estimator

    with plt.style.context(args.mplstylefile):
        fig, ax = plt.subplots(dpi=args.dpi)
        fig2, ax2 = plt.subplots(dpi=args.dpi)
        fig3, ax3 = plt.subplots(dpi=args.dpi)
        fig4, ax4 = plt.subplots(dpi=args.dpi)

        if args.binning_analysis:
            bfig, bax = plt.subplots(dpi=args.dpi)
            bfig2, bax2 = plt.subplots(dpi=args.dpi)
            bfig3, bax3 = plt.subplots(dpi=args.dpi)
            bfig4, bax4 = plt.subplots(dpi=args.dpi)


        for i,f in enumerate(args.filename):
            if args.labels:
                label = labels[i]
            else:
                label = os.path.basename(f)

            data = np.genfromtxt(f,names=custom_header,skip_header=args.skip_header,dtype=args.dtype,deletechars=args.deletechars)
            if args.estimator == "":
                msg = "Specify an estimator (--estimator {})".format(data.dtype.names)
                raise ValueError(msg)
            if not args.estimator in data.dtype.names:
                msg = "Estimator {} not in column headers {} for file {}".format(args.estimator, data.dtype.names,f)
                raise ValueError(msg)

            _avg, _err, fig, ax = averageplot_method(
                data[args.estimator],
                skip=args.skip,
                bin_size=args.bin_size,
                method=args.average_method,
                bootstrap_size=args.bootstrap_size,
                running_window=args.running_window,
                fig_ax=[fig,ax],
                label=label
            )

            fluctuations_avg, fluctuations_err, fig2, ax2 = averageplot_on_function_method(
                    fluctuations,data[args.estimator],data[args.estimator]**2,
                    skip=args.skip,
                    bin_size=args.bin_size,
                    method=args.average_method,
                    bootstrap_size=args.bootstrap_size,
                    running_window=args.running_window,
                    propogated_error_func=fluctuations_error,
                    fig_ax=[fig2,ax2],
                    label=label
            )

            skew_avg, skew_err, fig3, ax3 = averageplot_on_function_method(
                    skew,data[args.estimator],data[args.estimator]**2,data[args.estimator]**3,
                    skip=args.skip,
                    bin_size=args.bin_size,
                    method=args.average_method,
                    bootstrap_size=args.bootstrap_size,
                    running_window=args.running_window,
                    propogated_error_func=skew_error,
                    fig_ax=[fig3,ax3],
                    label=label
            )

            kurtosis_avg, kurtosis_err, fig4, ax4 = averageplot_on_function_method(
                    kurtosis,data[args.estimator],data[args.estimator]**2,data[args.estimator]**3,data[args.estimator]**4,
                    skip=args.skip,
                    bin_size=args.bin_size,
                    method=args.average_method,
                    bootstrap_size=args.bootstrap_size,
                    running_window=args.running_window,
                    propogated_error_func=kurtosis_error,
                    fig_ax=[fig4,ax4],
                    label=label
            )

            if args.binning_analysis:
                print("")
                print("")
                print(os.path.basename(f))
                print("Average")
                bfig, bax, estimator_autocorrelation_time, estimator_convergence_factor, estimator_max_bin_level, estimator_binned_avg, estimator_binned_err = binning_analysis(
                    data[args.estimator],
                    skip=args.skip,
                    method=args.average_method,
                    bootstrap_size=args.bootstrap_size,
                    running_window=args.running_window,
                    fig_ax=[bfig,bax],
                    label=label
                )
                print("")
                print("Fluctuations")

                bfig2, bax2, fluctuations_autocorrelation_time, fluctuations_convergence_factor, fluctuations_max_bin_level, fluctuations_binned_avg, fluctuations_binned_err = binning_analysis_on_function(
                    fluctuations, data[args.estimator], data[args.estimator]**2,
                    skip=args.skip,
                    method=args.average_method,
                    bootstrap_size=args.bootstrap_size,
                    running_window=args.running_window,
                    propogated_error_func=fluctuations_error,
                    fig_ax=[bfig2,bax2],
                    label=label
                )

                print("")
                print("Skew")
                bfig3, bax3, skew_autocorrelation_time, skew_convergence_factor, skew_max_bin_level, skew_binned_avg, skew_binned_err = binning_analysis_on_function(
                    skew, data[args.estimator], data[args.estimator]**2, data[args.estimator]**3,
                    skip=args.skip,
                    method=args.average_method,
                    bootstrap_size=args.bootstrap_size,
                    running_window=args.running_window,
                    propogated_error_func=skew_error,
                    fig_ax=[bfig3,bax3],
                    label=label
                )

                print("")
                print("Kurtosis")
                bfig4, bax4, kurtosis_autocorrelation_time, kurtosis_convergence_factor, kurtosis_max_bin_level, kurtosis_binned_avg, kurtosis_binned_err = binning_analysis_on_function(
                    kurtosis, data[args.estimator], data[args.estimator]**2, data[args.estimator]**3, data[args.estimator]**4,
                    skip=args.skip,
                    method=args.average_method,
                    bootstrap_size=args.bootstrap_size,
                    running_window=args.running_window,
                    propogated_error_func=kurtosis_error,
                    fig_ax=[bfig4,bax4],
                    label=label
                )

        ax.set_ylabel(r"$\langle {}\rangle$".format(pname))
        ax2.set_ylabel(r"$\langle {}^2\rangle - \langle {}\rangle^2$".format(pname,pname))
        ax3.set_ylabel(r"$\tilde{{\mu}}_{{3}}^{{{}}}$".format(pname))
        ax4.set_ylabel(r"$\kappa_{{{}}}$".format(pname))

        ax.set_title(r"Average")
        ax2.set_title(r"Fluctuations")
        ax3.set_title(r"Skew")
        ax4.set_title(r"Kurtosis")

        if args.legend:
            ax.legend()
            ax2.legend()
            ax3.legend()
            ax4.legend()

        if args.savefig:
            fig_fn = ("_{}_average".format(args.average_method)).join(os.path.splitext(args.savefig))
            fig2_fn = ("_{}_fluctuations".format(args.average_method)).join(os.path.splitext(args.savefig))
            fig3_fn = ("_{}_skew".format(args.average_method)).join(os.path.splitext(args.savefig))
            fig4_fn = ("_{}_kurtosis".format(args.average_method)).join(os.path.splitext(args.savefig))
            fig.savefig(fig_fn)
            fig2.savefig(fig2_fn)
            fig3.savefig(fig3_fn)
            fig4.savefig(fig4_fn)

        if args.binning_analysis:
            bax.set_ylabel(r"$\sigma_{{\langle {}\rangle}}$".format(pname))
            bax2.set_ylabel(r"$\sigma_{{\langle {}^2\rangle - \langle {}\rangle^2}}$".format(pname,pname))
            bax3.set_ylabel(r"$\sigma_{{\tilde{{\mu}}_{{3}}^{{{}}}}}$".format(pname))
            bax4.set_ylabel(r"$\sigma_{{\kappa_{{{}}}}}$".format(pname))

            try:
                tau_fmt = int(estimator_autocorrelation_time) + 1
            except:
                tau_fmt = "\mathrm{NaN}"
            try:
                test_nan = int(estimator_convergence_factor)
                cf_fmt = estimator_convergence_factor
            except:
                cf_fmt = "\mathrm{NaN}"
            _title = r"Average $\tau = {}$ $\mu_\mathrm{{C.F.}} = {}$".format(tau_fmt,cf_fmt)
            bax.set_title(_title)
            try:
                tau_fmt = int(fluctuations_autocorrelation_time) + 1
            except:
                tau_fmt = "\mathrm{NaN}"
            try:
                test_nan = int(fluctuations_convergence_factor)
                cf_fmt = fluctuations_convergence_factor
            except:
                cf_fmt = "\mathrm{NaN}"
            _title = r"Fluctuations $\tau = {}$ $\mu_\mathrm{{C.F.}} = {}$".format(tau_fmt,cf_fmt)
            bax2.set_title(_title)
            try:
                tau_fmt = int(skew_autocorrelation_time) + 1
            except:
                tau_fmt = "\mathrm{NaN}"
            try:
                test_nan = int(skew_convergence_factor)
                cf_fmt = skew_convergence_factor
            except:
                cf_fmt = "\mathrm{NaN}"
            _title = r"Skew $\tau = {}$ $\mu_\mathrm{{C.F.}} = {}$".format(tau_fmt,cf_fmt)
            bax3.set_title(_title)
            try:
                tau_fmt = int(kurtosis_autocorrelation_time) + 1
            except:
                tau_fmt = "\mathrm{NaN}"
            try:
                test_nan = int(kurtosis_convergence_factor)
                cf_fmt = kurtosis_convergence_factor
            except:
                cf_fmt = "\mathrm{NaN}"
            _title = r"Kurtosis $\tau = {}$ $\mu_\mathrm{{C.F.}} = {}$".format(tau_fmt,cf_fmt)
            bax4.set_title(_title)

            if args.legend:
                bax.legend()
                bax2.legend()
                bax3.legend()
                bax4.legend()

            if args.savefig:
                bfig_fn = ("_binning_{}_average".format(args.average_method)).join(os.path.splitext(args.savefig))
                bfig2_fn = ("_binning_{}_fluctuations".format(args.average_method)).join(os.path.splitext(args.savefig))
                bfig3_fn = ("_binning_{}_skew".format(args.average_method)).join(os.path.splitext(args.savefig))
                bfig4_fn = ("_binning_{}_kurtosis".format(args.average_method)).join(os.path.splitext(args.savefig))
                bfig.savefig(bfig_fn)
                bfig2.savefig(bfig2_fn)
                bfig3.savefig(bfig3_fn)
                bfig4.savefig(bfig4_fn)

    return 0



if __name__ == '__main__':
    sys.exit(main(sys.argv))

