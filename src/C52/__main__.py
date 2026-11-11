#!/usr/bin/env python3
"""module entry point"""

from argparse import ArgumentParser, Namespace
from sys import exit, stderr

from .primality import FermatTest, SolovayStrassenTest, MillerRabinTest, TestResultPrinter


OPTIMAL_BENCHMARK_ITER_AMT = 1000


def parse_args() -> Namespace:
    parser = ArgumentParser(
        prog='C52',
        description='CLI tool for probabilistic primality tests',
    )

    parser.add_argument('number', type=int)

    parser.add_argument(
        '-b', '--benchmark', dest='enable_benchmark',
        help='measure algorithm execution time',
        action='store_true',
    )
    parser.add_argument(
        '-s', '--stats', dest='display_stats',
        help='measure amt of positive results in 1000 runs',
        action='store_true',
    )
    parser.add_argument(
        '-v', '--verbose', dest='is_verbose',
        help='increase output verbosity',
        action='store_true',
    )

    parser.add_argument(
        '-r', '--rounds', dest='rounds',
        help='amount of algorithm rounds', type=int,
    )

    return parser.parse_args()


def main() -> None:
    try:
        args = parse_args()

        init_kwargs = {"rounds": args.rounds} if args.rounds is not None else {}
        fermat = FermatTest(**init_kwargs)
        solovay_strassen = SolovayStrassenTest(**init_kwargs)
        miller_rabin = MillerRabinTest(**init_kwargs)

        printer = TestResultPrinter(args.is_verbose)

        printer.print(fermat.is_prime(args.number))
        printer.print(solovay_strassen.is_prime(args.number))
        printer.print(miller_rabin.is_prime(args.number))

        if args.enable_benchmark:
            print('Benchmark:')
            print('    Fermat:', fermat.benchmark_ms(args.number, OPTIMAL_BENCHMARK_ITER_AMT) / OPTIMAL_BENCHMARK_ITER_AMT)
            print('    Solovay--Strassen:', solovay_strassen.benchmark_ms(args.number, OPTIMAL_BENCHMARK_ITER_AMT) / OPTIMAL_BENCHMARK_ITER_AMT)
            print('    Miller--Rabin:', miller_rabin.benchmark_ms(args.number, OPTIMAL_BENCHMARK_ITER_AMT) / OPTIMAL_BENCHMARK_ITER_AMT)

        if args.display_stats:
            print('Stats:')
            print('    Fermat:', fermat.gather_stats(args.number, OPTIMAL_BENCHMARK_ITER_AMT))
            print('    Solovay--Strassen:', solovay_strassen.gather_stats(args.number, OPTIMAL_BENCHMARK_ITER_AMT))
            print('    Miller--Rabin:', miller_rabin.gather_stats(args.number, OPTIMAL_BENCHMARK_ITER_AMT))

    except KeyboardInterrupt:
        print('\nCaught keyboard interrupt. Exiting...', file=stderr)
        exit(130)
    except Exception as e:
        print(f"E: {e}", file=stderr)
        exit(1)


if __name__ == '__main__':
    main()
