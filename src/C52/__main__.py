#!/usr/bin/env python3
"""module entry point"""

from argparse import ArgumentParser, Namespace
from sys import exit, stderr

from .primality import (
    FermatTest,
    SolovayStrassenTest,
    MillerRabinTest,
    TestResultPrinter,
    CarmichaelProvider,
)


OPTIMAL_BENCHMARK_ITER_AMT = 1000
DEFAULT_ROUND_AMT = 5


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
        '-c', '--carmichael', dest='carmichael_amt',
        help='run tests on carmichael numbers',
        type=int, default=0,
    )
    parser.add_argument(
        '-r', '--rounds', dest='rounds', default=DEFAULT_ROUND_AMT,
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

            fermat_ms = fermat.benchmark_ms(args.number, OPTIMAL_BENCHMARK_ITER_AMT)
            solovay_strassen_ms = solovay_strassen.benchmark_ms(args.number, OPTIMAL_BENCHMARK_ITER_AMT)
            miller_rabin_ms = miller_rabin.benchmark_ms(args.number, OPTIMAL_BENCHMARK_ITER_AMT)
            
            print('    Fermat:', fermat_ms / OPTIMAL_BENCHMARK_ITER_AMT)
            print('    Solovay--Strassen:', solovay_strassen_ms / OPTIMAL_BENCHMARK_ITER_AMT)
            print('    Miller--Rabin:', miller_rabin_ms / OPTIMAL_BENCHMARK_ITER_AMT)

        if args.display_stats:
            print('Stats:')
            print('    Fermat:', fermat.gather_stats(args.number, OPTIMAL_BENCHMARK_ITER_AMT))
            print('    Solovay--Strassen:', solovay_strassen.gather_stats(args.number, OPTIMAL_BENCHMARK_ITER_AMT))
            print('    Miller--Rabin:', miller_rabin.gather_stats(args.number, OPTIMAL_BENCHMARK_ITER_AMT))

        if args.carmichael_amt > 0:
            carmichael_provider = CarmichaelProvider()
            carmichael_nums = (carmichael_provider.get_number() for _ in range(args.carmichael_amt))

            for i in carmichael_nums:
                print('Carmichael number:', i.number)
                print('    Divisors:', *i.divisors)
                for j in (fermat, solovay_strassen, miller_rabin):
                    printer.print(j.is_prime(i.number))

    except KeyboardInterrupt:
        print('\nCaught keyboard interrupt. Exiting...', file=stderr)
        exit(130)
    except Exception as e:
        print(f"E: {e}", file=stderr)
        exit(1)


if __name__ == '__main__':
    main()
