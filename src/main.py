import argparse
from exons_coverage import CoverageManager
from coverage_seed import CoverageSeedManager
from cluster import ManagerCluster

class Manager:
    def __init__(self) -> None:
        pass

    def coverage(self, args):
        CoverageManager.coverage(args)

    def coverage_seed(self, args):
        CoverageSeedManager.coverage_seed(args)
        
    def cluster(self, args):
        ManagerCluster.cluster(args)
       

def main():
    manager = Manager()

    parser = argparse.ArgumentParser(description='Foo Bar')
    subparsers = parser.add_subparsers(dest='command', help='Commands to run', required=True)

    coverage = subparsers.add_parser('coverage', help='Получение покрытия для списка экзонов')
    coverage.add_argument('-o', '--output', help='Название для итогового файла с покрытием')
    coverage.add_argument("-b", "--bamfile", help="Bam file to get coverage", required=True)
    coverage.set_defaults(func=manager.coverage)

    coverage = subparsers.add_parser('cluster', help='Кластеризация покрытий экзонов относительно LINE')
    coverage.add_argument('-i', '--input', help='Название для итогового файла с покрытием', required=True)
    coverage.add_argument('-o', '--output', help='Название для итогового файла с кластеризованным покрытием', required=True)
    coverage.set_defaults(func=manager.cluster)

    coverage_seed = subparsers.add_parser('coverageseed', help='Получение покрытия затравок LINE')
    coverage_seed.add_argument('-i', '--input', help='Название для итогового файла с покрытием', required=True)
    coverage_seed.add_argument('-o', '--output', help='Название для итогового файла с кластеризованным покрытием', required=True)
    coverage_seed.add_argument("-b", "--bamfile", help="Bam file to get coverage", required=True)
    coverage_seed.set_defaults(func=manager.coverage_seed)



    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
  main()