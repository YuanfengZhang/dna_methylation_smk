from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from subprocess import run
from typing import Dict, List, Literal, Protocol
import polars as pl
from sample_sheet_parser import pipeline_decode, rev_pipeline_coding


@dataclass
class RepConfig:
    BaseName: str
    InDir: Path
    OutDir: Path


@dataclass
class SampleConfig:
    Counter: str
    Code: str
    Rep1: RepConfig
    Rep2: RepConfig


class ConvertCompressor(Protocol):
    def __call__(self, choice: Literal['A', 'B']) -> None:
        """
        定义一个函数协议，要求函数接受一个参数 choice
        其值必须是 'A' 或 'B'，并且返回值为 None
        """
        ...


class Convert2LazyFrame(Protocol):
    def __call__(self,
                 rep_config: RepConfig,
                 beta_col_name: str,
                 min_depth: int) -> pl.LazyFrame:
        ...


class BedGraphConverter(ABC):
    def __init__(self):
        """Abstract base class for BedGraph converters."""
        self._validate_subclass()

    @property
    @abstractmethod
    def Aim(self) -> str:
        """Aim of the conversion."""
        ...

    @property
    @abstractmethod
    def A(self) -> SampleConfig:
        """Sample A configuration."""
        ...

    @property
    @abstractmethod
    def B(self) -> SampleConfig:
        """Sample A configuration."""
        ...

    @property
    @abstractmethod
    def threads(self) -> int:
        """Number of threads to use for conversion."""
        ...

    def _validate_subclass(self):
        """Validate required attributes in subclasses."""
        required = ('Aim', 'A', 'B', 'threads')
        for attr in required:
            if not hasattr(self, attr):
                raise AttributeError(
                    f"Subclasses must define '{attr}' attribute"
                )

    @abstractmethod
    def convert(self) -> None:
        """Main conversion method to be implemented by subclasses."""
        ...

    @staticmethod
    def _sample_config(rep1_name: str, pipeline_code: str,
                       input_root_dir: Path, output_root_dir: Path) -> SampleConfig:

        rep2_name = rep1_name[:6] + '2' + rep1_name[7:]

        counter_dir = pipeline_decode(pipeline_code)
        counter = rev_pipeline_coding['Counter'][pipeline_code[-1]]

        input_dir1 = (input_root_dir / rep1_name / counter_dir)
        input_dir2 = (input_root_dir / rep2_name / counter_dir)

        if counter == 'dmrcate':
            out_dir1 = output_root_dir / 'dmrcate' / rep1_name
            out_dir2 = out_dir1
        else:
            out_dir1 = output_root_dir / rep1_name / counter_dir
            out_dir2 = output_root_dir / rep2_name / counter_dir

        return SampleConfig(Counter=counter,
                            Code=pipeline_code,
                            Rep1=RepConfig(BaseName=rep1_name, InDir=input_dir1, OutDir=out_dir1),
                            Rep2=RepConfig(BaseName=rep2_name, InDir=input_dir2, OutDir=out_dir2)
                            )


class BedGraph2Methylkit(BedGraphConverter):
    def __init__(self, s1r1_name: str, s2r1_name: str,
                 code1: str, code2: str,
                 input_root_dir: Path,
                 output_root_dir: Path,
                 threads: int = 1):
        self._aim: str = "methylkit"
        self._a: SampleConfig = self._sample_config(rep1_name=s1r1_name, pipeline_code=code1,
                                                    input_root_dir=input_root_dir,
                                                    output_root_dir=output_root_dir)
        self._b: SampleConfig = self._sample_config(rep1_name=s2r1_name, pipeline_code=code2,
                                                    input_root_dir=input_root_dir,
                                                    output_root_dir=output_root_dir)
        self._threads: int = max(threads, 1)
        super().__init__()

    @property
    def Aim(self) -> str:
        return self._aim

    @property
    def A(self) -> SampleConfig:
        return self._a

    @property
    def B(self) -> SampleConfig:
        return self._b

    @property
    def threads(self) -> int:
        return self._threads

    def convert(self):
        convert_options: Dict[str, ConvertCompressor] = {
            'astair': self._from_astair,
            'biscuit': self._from_biscuit,
            'bismark-bowtie2': self._from_bismark,
            'bismark-hisat2': self._from_bismark,
            'bsgenova': self._from_bsgenova,
            'fame': self._from_fame,
            'methyldackel': self._from_methyldackel
        }

        convert_options[self.A.Counter]('A')
        convert_options[self.A.Counter]('B')

    def _from_astair(self, choice: Literal['A', 'B']):
        if choice == 'A':
            sample = self.A
        else:
            sample = self.B

        r1_in: Path = sample.Rep1.InDir / f'{sample.Rep1.BaseName}.astair.mods.gz'
        r2_in: Path = sample.Rep2.InDir / f'{sample.Rep2.BaseName}.astair.mods.gz'
        r1_out: Path = sample.Rep1.OutDir / f'{sample.Rep1.BaseName}.{self.Aim}'
        r2_out: Path = sample.Rep2.OutDir / f'{sample.Rep2.BaseName}.{self.Aim}'

        for in_f, out_f in [(r1_in, r1_out), (r2_in, r2_out)]:
            if not in_f.exists():
                raise FileNotFoundError(in_f)

            (pl.scan_csv(in_f, separator='\t', has_header=True,
                         schema={'CHROM': pl.String, 'START': pl.Int64,
                                 'END': pl.Int64, 'MOD_LEVEL': pl.Float64,
                                 'MOD': pl.Int64, 'UNMOD': pl.Int64,
                                 'REF': pl.String, 'ALT': pl.String,
                                 'SPECIFIC_CONTEXT': pl.String, 'CONTEXT': pl.String,
                                 'SNV': pl.String, 'TOTAL_DEPTH': pl.Int64})
               .filter(pl.col('TOTAL_DEPTH') >= 5)
               .with_columns([pl.concat_str(pl.col('CHROM'),
                                            pl.col('START'),
                                            separator='.')
                                .alias('chrBase'),
                              pl.lit('F').alias('strand'),
                              (100 * pl.col('MOD') / pl.sum_horizontal('MOD', 'UNMOD'))
                             .alias('freqC'),
                              (100 * pl.col('UNMOD') / pl.sum_horizontal('MOD', 'UNMOD'))
                             .alias('freqT')])
               .select(['chrBase', 'CHROM', 'START', 'strand', 'TOTAL_DEPTH', 'freqC', 'freqT'])
               .rename({'CHROM': 'chr', 'START': 'base', 'TOTAL_DEPTH': 'coverage'})
               .sort(by=['chr', 'base'])
               .sink_csv(path=out_f, include_header=False, separator='\t',
                         sync_on_close='all', lazy=True)
               .collect())

            run(f'pigz -9 -p {self.threads} {out_f}', shell=True)

    def _from_biscuit(self, choice: Literal['A', 'B']):
        if choice == 'A':
            sample = self.A
        else:
            sample = self.B

        r1_in: Path = sample.Rep1.InDir / f'{sample.Rep1.BaseName}.epibed.gz'
        r2_in: Path = sample.Rep2.InDir / f'{sample.Rep2.BaseName}.epibed.gz'
        r1_out: Path = sample.Rep1.OutDir / f'{sample.Rep1.BaseName}.{self.Aim}'
        r2_out: Path = sample.Rep2.OutDir / f'{sample.Rep2.BaseName}.{self.Aim}'

        for in_f, out_f in [(r1_in, r1_out), (r2_in, r2_out)]:
            if not in_f.exists():
                raise FileNotFoundError(in_f)

            (pl.scan_csv(in_f, separator='\t', has_header=False,
                         schema={'CHROM': pl.String, 'START': pl.Int64,
                                 'END': pl.Int64, 'RefBase': pl.String,
                                 'CytosineGroup': pl.String, '2Base': pl.String, '5Base': pl.String,
                                 'freqC': pl.Float64, 'coverage': pl.Int64})
               .filter(pl.col('coverage') >= 5)
               .with_columns(100 * pl.col('freqC'))
               .with_columns([pl.concat_str(pl.col('CHROM'),
                                            pl.col('START'),
                                            separator='.')
                                .alias('chrBase'),
                              pl.lit('F').alias('strand'),
                              (100 - pl.col('freqC')).alias('freqT')])
               .select(['chrBase', 'CHROM', 'START', 'strand', 'coverage', 'freqC', 'freqT'])
               .rename({'CHROM': 'chr', 'START': 'base'})
               .sort(by=['chr', 'base'])
               .sink_csv(path=out_f, include_header=False, separator='\t',
                         sync_on_close='all', lazy=True)
               .collect())

            run(f'pigz -9 -p {self.threads} {out_f}', shell=True)

    def _from_bismark(self, choice: Literal['A', 'B']):
        # TODO: FINISH THIS
        pass

    def _from_bsgenova(self, choice: Literal['A', 'B']):
        if choice == 'A':
            sample = self.A
        else:
            sample = self.B

        r1_in: Path = sample.Rep1.InDir / f'{sample.Rep1.BaseName}.ATCGmap.gz'
        r2_in: Path = sample.Rep2.InDir / f'{sample.Rep2.BaseName}.ATCGmap.gz'
        r1_out: Path = sample.Rep1.OutDir / f'{sample.Rep1.BaseName}.{self.Aim}'
        r2_out: Path = sample.Rep2.OutDir / f'{sample.Rep2.BaseName}.{self.Aim}'

        for in_f, out_f in [(r1_in, r1_out), (r2_in, r2_out)]:
            if not in_f.exists():
                raise FileNotFoundError(in_f)

            (pl.scan_csv(in_f, separator='\t', has_header=True)
               .filter(pl.col('depth') >= 5)
               .with_columns(100 * pl.col('beta').alias('freqC'))
               .with_columns([pl.concat_str(pl.col('chr'), pl.col('pos'),
                                            separator='.')
                                .alias('chrBase'),
                              (100 - pl.col('freqC')).alias('freqT')])
               .with_columns(pl.when(pl.col('base') == 'C')
                               .then(pl.lit('F'))
                               .otherwise(pl.lit('R'))
                               .alias('strand'))
               .select(['chrBase', 'chr', 'pos', 'strand', 'depth', 'freqC', 'freqT'])
               .rename({'pos': 'base'})
               .sort(by=['chr', 'base'])
               .sink_csv(path=out_f, include_header=False, separator='\t',
                         sync_on_close='all', lazy=True)
               .collect())
            run(f'pigz -9 -p {self.threads} {out_f}', shell=True)

    def _from_fame(self, choice: Literal['A', 'B']):
        if choice == 'A':
            sample = self.A
        else:
            sample = self.B

        r1_in: Path = sample.Rep1.InDir / f'{sample.Rep1.BaseName}.bedgraph.zst'
        r2_in: Path = sample.Rep2.InDir / f'{sample.Rep2.BaseName}.bedgraph.zst'
        r1_out: Path = sample.Rep1.OutDir / f'{sample.Rep1.BaseName}.{self.Aim}'
        r2_out: Path = sample.Rep2.OutDir / f'{sample.Rep2.BaseName}.{self.Aim}'

        for in_f, out_f in [(r1_in, r1_out), (r2_in, r2_out)]:
            if not in_f.exists():
                raise FileNotFoundError(in_f)

            (pl.scan_csv(in_f, separator='\t', has_header=True)
               .rename({'start': 'base', 'beta': 'freqC'})
               .with_columns(pl.sum_horizontal('m+', 'm-', 'u+', 'u-').alias('depth'))
               .filter(pl.col('depth') >= 5)
               .with_columns([pl.concat_str(pl.col('chrom'), pl.col('base'),
                                            separator='.')
                                .alias('chrBase'),
                              (100 - pl.col('freqC')).alias('freqT'),
                              pl.lit('F').alias('strand')])
               .select(['chrBase', 'chrom', 'base', 'strand', 'depth', 'freqC', 'freqT'])
               .sort(by=['chrom', 'base'])
               .sink_csv(path=out_f, include_header=False, separator='\t',
                         sync_on_close='all', lazy=True)
               .collect())
            run(f'pigz -9 -p {self.threads} {out_f}', shell=True)

    def _from_methyldackel(self, choice: Literal['A', 'B']):
        if choice == 'A':
            sample = self.A
        else:
            sample = self.B

        r1_in: Path = sample.Rep1.InDir / f'{sample.Rep1.BaseName}.merged.bedgraph.gz'
        r2_in: Path = sample.Rep2.InDir / f'{sample.Rep2.BaseName}.merged.bedgraph.gz'
        r1_out: Path = sample.Rep1.OutDir / f'{sample.Rep1.BaseName}.{self.Aim}'
        r2_out: Path = sample.Rep2.OutDir / f'{sample.Rep2.BaseName}.{self.Aim}'

        for in_f, out_f in [(r1_in, r1_out), (r2_in, r2_out)]:
            if not in_f.exists():
                raise FileNotFoundError(in_f)

            (pl.scan_csv(in_f, separator='\t', has_header=False,
                         schema={'chrom': pl.String, 'base': pl.Int64,
                                 'end': pl.Int64, 'freqC': pl.Float64,
                                 'm': pl.Int64, 'u': pl.Int64})
             .with_columns(pl.sum_horizontal('m', 'u').alias('depth'))
             .filter(pl.col('depth') >= 5)
             .with_columns([pl.concat_str(pl.col('chrom'), pl.col('base'),
                                          separator='.')
                           .alias('chrBase'),
                            (100 - pl.col('freqC')).alias('freqT'),
                            pl.lit('F').alias('strand')])
             .select(['chrBase', 'chrom', 'base', 'strand', 'depth', 'freqC', 'freqT'])
             .sort(by=['chrom', 'base'])
             .sink_csv(path=out_f, include_header=False, separator='\t',
                       sync_on_close='all', lazy=True)
             .collect())
            run(f'pigz -9 -p {self.threads} {out_f}', shell=True)


class BedGraph2Methylsig(BedGraphConverter):
    def __init__(self, s1r1_name: str, s2r1_name: str,
                 code1: str, code2: str,
                 input_root_dir: Path,
                 output_root_dir: Path,
                 threads: int = 1):
        self._aim: str = "methylkit"
        self._a: SampleConfig = self._sample_config(rep1_name=s1r1_name, pipeline_code=code1,
                                                    input_root_dir=input_root_dir,
                                                    output_root_dir=output_root_dir)
        self._b: SampleConfig = self._sample_config(rep1_name=s2r1_name, pipeline_code=code2,
                                                    input_root_dir=input_root_dir,
                                                    output_root_dir=output_root_dir)
        self._threads: int = max(threads, 1)
        super().__init__()

    @property
    def Aim(self) -> str:
        return self._aim

    @property
    def A(self) -> SampleConfig:
        return self._a

    @property
    def B(self) -> SampleConfig:
        return self._b

    @property
    def threads(self) -> int:
        return self._threads

    def convert(self):
        convert_options: Dict[str, ConvertCompressor] = {
            'astair': self._from_astair,
            'biscuit': self._from_biscuit,
            'bismark-bowtie2': self._from_bismark,
            'bismark-hisat2': self._from_bismark,
            'bsgenova': self._from_bsgenova,
            'fame': self._from_fame,
            'methyldackel': self._from_methyldackel
        }

        convert_options[self.A.Counter]('A')
        convert_options[self.A.Counter]('B')

    def _from_astair(self, choice: Literal['A', 'B']):
        if choice == 'A':
            sample = self.A
        else:
            sample = self.B

        r1_in: Path = sample.Rep1.InDir / f'{sample.Rep1.BaseName}.astair.mods.gz'
        r2_in: Path = sample.Rep2.InDir / f'{sample.Rep2.BaseName}.astair.mods.gz'
        r1_out: Path = sample.Rep1.OutDir / f'{sample.Rep1.BaseName}.{self.Aim}'
        r2_out: Path = sample.Rep2.OutDir / f'{sample.Rep2.BaseName}.{self.Aim}'

        for in_f, out_f in [(r1_in, r1_out), (r2_in, r2_out)]:
            if not in_f.exists():
                raise FileNotFoundError(in_f)

            (pl.scan_csv(in_f, separator='\t', has_header=True,
                         schema={'CHROM': pl.String, 'START': pl.Int64,
                                 'END': pl.Int64, 'MOD_LEVEL': pl.Float64,
                                 'MOD': pl.Int64, 'UNMOD': pl.Int64,
                                 'REF': pl.String, 'ALT': pl.String,
                                 'SPECIFIC_CONTEXT': pl.String, 'CONTEXT': pl.String,
                                 'SNV': pl.String, 'TOTAL_DEPTH': pl.Int64})
             .filter(pl.col('TOTAL_DEPTH') >= 5)
             .with_columns(pl.col('MOD_LEVEL') * 100)
             .select(['CHROM', 'START', 'MOD_LEVEL', 'TOTAL_DEPTH'])
             .sort(by=['CHROM', 'START'])
             .sink_csv(path=out_f, include_header=False, separator='\t',
                       sync_on_close='all', lazy=True)
             .collect())

            run(f'zstd -19 --rm -T{self.threads} {out_f}', shell=True)

    def _from_biscuit(self, choice: Literal['A', 'B']):
        if choice == 'A':
            sample = self.A
        else:
            sample = self.B

        r1_in: Path = sample.Rep1.InDir / f'{sample.Rep1.BaseName}.epibed.gz'
        r2_in: Path = sample.Rep2.InDir / f'{sample.Rep2.BaseName}.epibed.gz'
        r1_out: Path = sample.Rep1.OutDir / f'{sample.Rep1.BaseName}.{self.Aim}'
        r2_out: Path = sample.Rep2.OutDir / f'{sample.Rep2.BaseName}.{self.Aim}'

        for in_f, out_f in [(r1_in, r1_out), (r2_in, r2_out)]:
            if not in_f.exists():
                raise FileNotFoundError(in_f)

            (pl.scan_csv(in_f, separator='\t', has_header=False,
                         schema={'CHROM': pl.String, 'START': pl.Int64,
                                 'END': pl.Int64, 'RefBase': pl.String,
                                 'CytosineGroup': pl.String, '2Base': pl.String, '5Base': pl.String,
                                 'freqC': pl.Float64, 'coverage': pl.Int64})
               .filter(pl.col('coverage') >= 5)
               .with_columns((pl.col('freqC') * 100).alias('beta'))
               .select(['CHROM', 'START', 'beta', 'coverage'])
               .sort(by=['CHROM', 'START'])
               .sink_csv(path=out_f, include_header=False, separator='\t',
                         sync_on_close='all', lazy=True)
               .collect())

            run(f'zstd -19 --rm -T{self.threads} {out_f}', shell=True)

    def _from_bismark(self, choice: Literal['A', 'B']):
        # TODO: FINISH THIS
        pass

    def _from_bsgenova(self, choice: Literal['A', 'B']):
        if choice == 'A':
            sample = self.A
        else:
            sample = self.B

        r1_in: Path = sample.Rep1.InDir / f'{sample.Rep1.BaseName}.ATCGmap.gz'
        r2_in: Path = sample.Rep2.InDir / f'{sample.Rep2.BaseName}.ATCGmap.gz'
        r1_out: Path = sample.Rep1.OutDir / f'{sample.Rep1.BaseName}.{self.Aim}'
        r2_out: Path = sample.Rep2.OutDir / f'{sample.Rep2.BaseName}.{self.Aim}'

        for in_f, out_f in [(r1_in, r1_out), (r2_in, r2_out)]:
            if not in_f.exists():
                raise FileNotFoundError(in_f)

            (pl.scan_csv(in_f, separator='\t', has_header=True)
               .filter(pl.col('depth') >= 5)
               .with_columns([(pl.col('pos') - 2).alias('start'),
                              100 * pl.col('beta').alias('freqC')])
               .select(['chr', 'start', 'freqC', 'depth'])
               .sort(by=['chr', 'start'])
               .sink_csv(path=out_f, include_header=False, separator='\t',
                         sync_on_close='all', lazy=True)
               .collect())
            run(f'zstd -19 --rm -T{self.threads} {out_f}', shell=True)

    def _from_fame(self, choice: Literal['A', 'B']):
        if choice == 'A':
            sample = self.A
        else:
            sample = self.B

        r1_in: Path = sample.Rep1.InDir / f'{sample.Rep1.BaseName}.bedgraph.zst'
        r2_in: Path = sample.Rep2.InDir / f'{sample.Rep2.BaseName}.bedgraph.zst'
        r1_out: Path = sample.Rep1.OutDir / f'{sample.Rep1.BaseName}.{self.Aim}'
        r2_out: Path = sample.Rep2.OutDir / f'{sample.Rep2.BaseName}.{self.Aim}'

        for in_f, out_f in [(r1_in, r1_out), (r2_in, r2_out)]:
            if not in_f.exists():
                raise FileNotFoundError(in_f)

            (pl.scan_csv(in_f, separator='\t', has_header=True)
               .with_columns([pl.sum_horizontal('m+', 'm-', 'u+', 'u-').alias('depth')])
               .filter(pl.col('depth') >= 5)
               .select(['chrom', 'start', 'beta', 'depth'])
               .sort(by=['chrom', 'start'])
               .sink_csv(path=out_f, include_header=False, separator='\t',
                         sync_on_close='all', lazy=True)
               .collect())
            run(f'zstd -19 --rm -T{self.threads} {out_f}', shell=True)

    def _from_methyldackel(self, choice: Literal['A', 'B']):
        if choice == 'A':
            sample = self.A
        else:
            sample = self.B

        r1_in: Path = sample.Rep1.InDir / f'{sample.Rep1.BaseName}.merged.bedgraph.gz'
        r2_in: Path = sample.Rep2.InDir / f'{sample.Rep2.BaseName}.merged.bedgraph.gz'
        r1_out: Path = sample.Rep1.OutDir / f'{sample.Rep1.BaseName}.{self.Aim}'
        r2_out: Path = sample.Rep2.OutDir / f'{sample.Rep2.BaseName}.{self.Aim}'

        for in_f, out_f in [(r1_in, r1_out), (r2_in, r2_out)]:
            if not in_f.exists():
                raise FileNotFoundError(in_f)

            (pl.scan_csv(in_f, separator='\t', has_header=False,
                         schema={'chrom': pl.String, 'start': pl.Int64,
                                 'end': pl.Int64, 'beta': pl.Float64,
                                 'm': pl.Int64, 'u': pl.Int64})
             .with_columns(pl.sum_horizontal('m', 'u').alias('depth'))
             .filter(pl.col('depth') >= 5)
             .select(['chrom', 'start', 'beta', 'depth'])
             .sort(by=['chrom', 'start'])
             .sink_csv(path=out_f, include_header=False, separator='\t',
                       sync_on_close='all', lazy=True)
             .collect())
            run(f'zstd -19 --rm -T{self.threads} {out_f}', shell=True)


class BedGraph2Dmrcate(BedGraph2Methylsig):
    def __init__(self, s1r1_name: str, s2r1_name: str, code1: str, code2: str,
                 input_root_dir: Path, output_root_dir: Path, threads: int = 1):
        super().__init__(s1r1_name=s1r1_name, s2r1_name=s2r1_name,
                         code1=code1, code2=code2, threads=threads,
                         input_root_dir=input_root_dir, output_root_dir=output_root_dir)
        self._aim = "dmrcate"


class BedGraph2Metilene(BedGraphConverter):
    def __init__(self, s1r1_name: str, s2r1_name: str,
                 code1: str, code2: str,
                 input_root_dir: Path,
                 output_root_dir: Path,
                 threads: int = 1):
        self._aim: str = "metilene"
        self._a: SampleConfig = self._sample_config(rep1_name=s1r1_name, pipeline_code=code1,
                                                    input_root_dir=input_root_dir,
                                                    output_root_dir=output_root_dir)
        self._b: SampleConfig = self._sample_config(rep1_name=s2r1_name, pipeline_code=code2,
                                                    input_root_dir=input_root_dir,
                                                    output_root_dir=output_root_dir)
        self._threads: int = max(threads, 1)
        super().__init__()

    @property
    def Aim(self) -> str:
        return self._aim

    @property
    def A(self) -> SampleConfig:
        return self._a

    @property
    def B(self) -> SampleConfig:
        return self._b

    @property
    def threads(self) -> int:
        return self._threads

    @staticmethod
    def flip_basename(basename: str) -> str:
        parts: List[str] = basename.split('_')
        return '_'.join((parts[0], parts[3], parts[1], parts[2]))

    def convert(self):
        convert_options: Dict[str, Convert2LazyFrame] = {
            'astair': self._from_astair,
            'biscuit': self._from_biscuit,
            'bismark-bowtie2': self._from_bismark,
            'bismark-hisat2': self._from_bismark,
            'bsgenova': self._from_bsgenova,
            'fame': self._from_fame,
            'methyldackel': self._from_methyldackel
        }

        (s1r1_l,
         s1r2_l,
         s2r1_l,
         s2r2_l) = [self.flip_basename(_n) for _n in [self.A.Rep1.BaseName, self.A.Rep2.BaseName,
                                                      self.B.Rep1.BaseName, self.B.Rep2.BaseName]]

        df_s1r1 = convert_options[self.A.Counter](rep_config=self.A.Rep1, beta_col_name=s1r1_l, min_depth=5)
        df_s1r2 = convert_options[self.A.Counter](rep_config=self.A.Rep2, beta_col_name=s1r2_l, min_depth=5)
        df_s2r1 = convert_options[self.B.Counter](rep_config=self.B.Rep1, beta_col_name=s2r1_l, min_depth=5)
        df_s2r2 = convert_options[self.B.Counter](rep_config=self.B.Rep2, beta_col_name=s2r2_l, min_depth=5)

        (df_s1r1.join(other=df_s1r2, on=['chrom', 'start'], how='inner'))

        solo_out: Path = (self.A.Rep1.OutDir /  # noqa
                          f'{self.A.Rep1.BaseName}_vs_{self.B.Rep1.BaseName}' /  # noqa
                          f'{self.A.Code}.{self.B.Code}.{self.Aim}')

        (df_s1r1.join(other=df_s1r2, on=['chrom', 'start'], how='inner')
                .join(other=df_s2r1, on=['chrom', 'start'], how='inner')
                .join(other=df_s2r2, on=['chrom', 'start'], how='inner')
                .drop_nulls()
                .sink_csv(path=solo_out, include_header=False, separator='\t',
                          sync_on_close='all', lazy=True))

    @staticmethod
    def _from_astair(rep_config: RepConfig, beta_col_name: str, min_depth: int = 5) -> pl.LazyFrame:
        in_f: Path = rep_config.InDir / f'{rep_config.BaseName}.astair.mods.gz'

        if not in_f.exists():
            raise FileNotFoundError(in_f)

        return (pl.scan_csv(in_f, separator='\t', has_header=True,
                            schema={'CHROM': pl.String, 'START': pl.Int64,
                                    'END': pl.Int64, 'MOD_LEVEL': pl.Float64,
                                    'MOD': pl.Int64, 'UNMOD': pl.Int64,
                                    'REF': pl.String, 'ALT': pl.String,
                                    'SPECIFIC_CONTEXT': pl.String, 'CONTEXT': pl.String,
                                    'SNV': pl.String, 'TOTAL_DEPTH': pl.Int64})
                  .filter(pl.col('TOTAL_DEPTH') >= min_depth)
                  .with_columns(pl.col('MOD_LEVEL') * 100)
                  .rename({'CHROM': 'chrom', 'START': 'start', 'MOD_LEVEL': beta_col_name})
                  .select(['chrom', 'start', beta_col_name])
                  .unique(subset=['chrom', 'start'])
                  .sort(by=['chrom', 'start']))

    @staticmethod
    def _from_biscuit(rep_config: RepConfig, beta_col_name: str, min_depth: int = 5) -> pl.LazyFrame:
        in_f: Path = rep_config.InDir / f'{rep_config.BaseName}.epibed.gz'

        if not in_f.exists():
            raise FileNotFoundError(in_f)

        return (pl.scan_csv(in_f, separator='\t', has_header=False,
                            schema={'CHROM': pl.String, 'START': pl.Int64,
                                    'END': pl.Int64, 'RefBase': pl.String,
                                    'CytosineGroup': pl.String, '2Base': pl.String, '5Base': pl.String,
                                    'freqC': pl.Float64, 'coverage': pl.Int64})
                  .filter(pl.col('coverage') >= min_depth)
                  .with_columns((pl.col('freqC') * 100).alias(beta_col_name))
                  .select(['CHROM', 'START', beta_col_name])
                  .rename({'CHROM': 'chrom', 'START': 'start'})
                  .unique(subset=['chrom', 'start'])
                  .sort(by=['chrom', 'start']))

    @staticmethod
    def _from_bismark(rep_config: RepConfig, beta_col_name: str, min_depth: int = 5) -> pl.LazyFrame:
        # TODO: FINISH THIS
        in_f: Path = rep_config.InDir / f'{rep_config.BaseName}.bismark.cov.gz'
        if not in_f.exists():
            raise FileNotFoundError(in_f)

        return (pl.scan_csv(in_f, separator='\t', has_header=True))

    @staticmethod
    def _from_bsgenova(rep_config: RepConfig, beta_col_name: str, min_depth: int = 5) -> pl.LazyFrame:
        in_f: Path = rep_config.InDir / f'{rep_config.BaseName}.ATCGmap.gz'

        if not in_f.exists():
            raise FileNotFoundError(in_f)

        return (pl.scan_csv(in_f, separator='\t', has_header=True)
                  .filter(pl.col('depth') >= min_depth)
                  .with_columns([(pl.col('pos') - 2).alias('start'),
                                 100 * pl.col('beta').alias(beta_col_name)])
                  .select(['chr', 'start', beta_col_name])
                  .rename({'chr': 'chrom'})
                  .unique(subset=['chrom', 'start'])
                  .sort(by=['chrom', 'start']))

    @staticmethod
    def _from_fame(rep_config: RepConfig, beta_col_name: str, min_depth: int = 5) -> pl.LazyFrame:
        in_f: Path = rep_config.InDir / f'{rep_config.BaseName}.ATCGmap.gz'

        if not in_f.exists():
            raise FileNotFoundError(in_f)

        return (pl.scan_csv(in_f, separator='\t', has_header=True)
                  .filter(pl.sum_horizontal('m+', 'm-', 'u+', 'u-') >= min_depth)
                  .rename({'beta': beta_col_name})
                  .select(['chrom', 'start', beta_col_name])
                  .unique(subset=['chrom', 'start'])
                  .sort(by=['chrom', 'start']))

    @staticmethod
    def _from_methyldackel(rep_config: RepConfig, beta_col_name: str, min_depth: int = 5) -> pl.LazyFrame:
        in_f: Path = rep_config.InDir / f'{rep_config.BaseName}.merged.bedgraph.gz'

        if not in_f.exists():
            raise FileNotFoundError(in_f)

        return (pl.scan_csv(in_f, separator='\t', has_header=False,
                            schema={'chrom': pl.String, 'start': pl.Int64,
                                    'end': pl.Int64, 'beta': pl.Float64,
                                    'm': pl.Int64, 'u': pl.Int64})
                  .filter(pl.sum_horizontal('m', 'u') >= min_depth)
                  .rename({'beta': beta_col_name})
                  .select(['chrom', 'start', beta_col_name])
                  .unique(subset=['chrom', 'start'])
                  .sort(by=['chrom', 'start']))


if __name__ == '__main__':
    # TODO: TEST ALL
    pass
