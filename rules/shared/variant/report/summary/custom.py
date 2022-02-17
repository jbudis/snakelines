import os
import re

from bs4 import BeautifulSoup
from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML


class Status:
    PASS = 'PASS'
    WARN = 'WARN'
    FAIL = 'FAIL'
    NONE = 'MISSING'
    
    STATUS_LIST = [PASS, WARN, FAIL, NONE]
    
    STATUS_SCORE = {
        NONE: 0,
        PASS: 1,
        WARN: 2,
        FAIL: 3
    }


class FastQCParser:
    PER_BASE_QUALITY = 'Per base sequence quality'
    PER_TILE_QUALITY = 'Per tile sequence quality'
    PER_SEQUENCE_QUALITY = 'Per sequence quality scores'
    PER_SEQUENCE_GC_CONTENT = 'Per sequence GC content'
    DUPLICATION_LEVELS = 'Sequence Duplication Levels'
    PER_BASE_N_CONTENT = 'Per base N content'
    ADAPTER_CONTENT = 'Adapter Content'
    
    VALIDATORS = [
        PER_BASE_QUALITY,
        PER_TILE_QUALITY,
        PER_SEQUENCE_QUALITY,
        PER_SEQUENCE_GC_CONTENT,
        DUPLICATION_LEVELS,
        PER_BASE_N_CONTENT,
        ADAPTER_CONTENT
    ]
    
    @classmethod
    def parse_summary(cls, filename: str) -> list:
        summary = []
        summary_validators = {}
        with open(filename, 'r') as summary_file:
            for line in summary_file:
                status, name, _ = line.split('\t')
                assert status in Status.STATUS_LIST
                summary_validators[name] = status
        
        for name in cls.VALIDATORS:
            if name in summary_validators:
                summary.append({'name': name, 'status': summary_validators[name]})
            else:
                summary.append({'name': name, 'status': Status.NONE})
        
        return summary
    
    @classmethod
    def join_summary_pair(cls, summary_r1: list, summary_r2: list) -> list:
        assert len(summary_r1) == len(summary_r2)
        
        summary = []
        for i in range(len(summary_r1)):
            status_r1 = summary_r1[i]['status']
            status_r2 = summary_r2[i]['status']
            
            assert status_r1 in Status.STATUS_LIST
            assert status_r2 in Status.STATUS_LIST
            
            # select worse status
            if Status.STATUS_SCORE[status_r1] >= Status.STATUS_SCORE[status_r2]:
                summary.append(summary_r1[i])
            else:
                summary.append(summary_r2[i])
        
        return summary
    
    @classmethod
    def parse_summary_pair(cls, original_fastqc_r1: str, original_fastqc_r2: str) -> list:
        summary_r1 = cls.parse_summary(original_fastqc_r1)
        summary_r2 = cls.parse_summary(original_fastqc_r2)
        return cls.join_summary_pair(summary_r1, summary_r2)


class QualimapParser:
    ALL_MAPPED_RATIO = 'All mapped reads'
    MAPPED_RATIO = 'Panel mapped reads'
    MAPPING_QUALITY = 'Mean mapping quality'
    ERROR_RATE = 'General error rate'
    
    MEAN_COVERAGE = 'Mean coverage'
    COVERAGE_RATIO = 'Covered reads ratio'
    MULTI_COVERAGE_RATIO = '10-fold covered reads ratio'
    
    FLOAT_REGEX = '([0-9]*\.)?[0-9]+'
    PERCENTAGE_REGEX = FLOAT_REGEX + '%'
    MEAN_COVERAGE_REGEX = 'mean coverageData = ' + FLOAT_REGEX + 'X'
    COVERAGE_RATIO_REGEX = 'There is a ' + PERCENTAGE_REGEX + ' of reference with a coverageData >='
    ERROR_RATE_REGEX = 'general error rate = ' + FLOAT_REGEX
    MAPPING_QUALITY_REGEX = 'mean mapping quality = ' + FLOAT_REGEX
    
    @classmethod
    def _parse_value_from_line(cls, filename, line_regex: str, is_percentage: bool = False) -> float:
        result = 0.0
        with open(filename, 'r') as in_file:
            for line in in_file:
                stripped_line = line.strip()
                if re.match(line_regex, stripped_line):
                    result = cls._parse_value(stripped_line, is_percentage)
                    break
        
        return result
    
    @classmethod
    def _parse_value(cls, text, is_percentage=False) -> float:
        """
        :param text:
        :param is_percentage: if text value is in form of percetnage (e.g. '98.3 %')
        :return: percentage is returned as ratio (e.g. '98.3 %' == 0.983)
        """
        if is_percentage:
            match_obj = re.search(cls.PERCENTAGE_REGEX, text)
        else:
            match_obj = re.search(cls.FLOAT_REGEX, text)
        
        if match_obj is not None:
            if is_percentage:
                result = float(match_obj.group()[:-1]) / 100
            else:
                result = float(match_obj.group())
        else:
            raise ValueError('percentage not found in text "%s"' % text)
        
        return result
    
    @classmethod
    def _coverage_ratio(cls, results_filename: str, from_coverage: int = 1) -> float:
        regex = cls.COVERAGE_RATIO_REGEX + ' ' + str(from_coverage) + 'X'
        result = cls._parse_value_from_line(results_filename, regex, True)
        return result
    
    @classmethod
    def _mean_coverage(cls, results_filename: str) -> float:
        result = cls._parse_value_from_line(results_filename, cls.MEAN_COVERAGE_REGEX)
        return result
    
    @classmethod
    def _parse_general_error(cls, results_filename: str) -> float:
        result = cls._parse_value_from_line(results_filename, cls.ERROR_RATE_REGEX)
        return result
    
    @classmethod
    def _parse_mapping_quality(cls, results_filename: str) -> float:
        result = cls._parse_value_from_line(results_filename, cls.MAPPING_QUALITY_REGEX)
        return result
    
    @classmethod
    def _all_mapped_reads(cls, report_filename: str) -> float:
        with open(report_filename, 'r') as in_file:
            html = in_file.read()
        
        soup = BeautifulSoup(html, 'html.parser')
        element = soup \
            .find('h3', text='Globals') \
            .find_next_sibling('table') \
            .find('td', text='Mapped reads') \
            .find_next_sibling()
        
        return cls._parse_value(element.text, True)
    
    @classmethod
    def _panel_mapped_reads(cls, report_filename: str) -> float:
        result = None
        
        with open(report_filename, 'r') as in_file:
            html = in_file.read()
        
        soup = BeautifulSoup(html, 'html.parser')
        element = soup.find('h3', text='Globals (inside of regions)')
        if element is not None:
            value_el = element.find_next_sibling('table') \
                .find('td', text='Mapped reads') \
                .find_next_sibling()
            result = cls._parse_value(value_el.text, True)
        
        return result
    
    @classmethod
    def parse_summmary(cls, qualimap_txt: str, report_html: str) -> list:
        summary = []
        
        all_mapped_ratio = cls._all_mapped_reads(report_html)
        status = cls._value2status(all_mapped_ratio, 0.9, 0.95)
        summary.append({
            'name': cls.ALL_MAPPED_RATIO,
            'status': status,
            'value': all_mapped_ratio,
            'value_str': '%.2f%%' % (all_mapped_ratio * 100)
        })
        
        mapped_ratio = cls._panel_mapped_reads(report_html)
        if mapped_ratio is not None:
            # panel was used
            status = cls._value2status(mapped_ratio, 0.6, 0.8)
            summary.append({
                'name': cls.MAPPED_RATIO,
                'status': status,
                'value': mapped_ratio,
                'value_str': '%.2f%%' % (mapped_ratio * 100)
            })
        
        mean_coverage = cls._mean_coverage(qualimap_txt)
        status = cls._value2status(mean_coverage, 0.2, 0.3)
        summary.append({
            'name': cls.MEAN_COVERAGE,
            'status': status,
            'value': mean_coverage,
            'value_str': '%.2fX' % mean_coverage,
        })
        
        mapping_quality = cls._parse_mapping_quality(qualimap_txt)
        status = cls._value2status(mapping_quality, 0.35, 0.38)
        summary.append({
            'name': cls.MAPPING_QUALITY,
            'status': status,
            'value': mapping_quality,
            'value_str': '%.2f' % mapping_quality,
        })
        
        general_error = cls._parse_general_error(qualimap_txt)
        status = cls._value2status(general_error, 0.5, 0.4)
        summary.append({
            'name': cls.ERROR_RATE,
            'status': status,
            'value': general_error,
            'value_str': '%.4f' % general_error,
        })
        
        covered_ratio = cls._coverage_ratio(qualimap_txt, 1)
        status = cls._value2status(covered_ratio, 0.8, 0.9)
        summary.append({
            'name': cls.COVERAGE_RATIO,
            'status': status,
            'value': covered_ratio,
            'value_str': '%.2f%%' % (covered_ratio * 100),
        })
        
        multi_covered_ratio = cls._coverage_ratio(qualimap_txt, 10)
        status = cls._value2status(multi_covered_ratio, 0.7, 0.8)
        summary.append({
            'name': cls.MULTI_COVERAGE_RATIO,
            'status': status,
            'value': multi_covered_ratio,
            'value_str': '%.2f%%' % (multi_covered_ratio * 100),
        })
        
        return summary
    
    @classmethod
    def _value2status(cls, value: float, fail_limit: float, warn_limit: float) -> str:
        assert fail_limit != warn_limit
        
        status = Status.FAIL
        if warn_limit > fail_limit:
            # greater is better
            if value >= warn_limit:
                status = Status.PASS
            elif value >= fail_limit:
                status = Status.WARN
        
        elif fail_limit > warn_limit:
            # lower is better
            if value <= warn_limit:
                status = Status.PASS
            elif value <= fail_limit:
                status = Status.WARN
        
        return status


class GatkCallingParser:
    HEADER_LINE_ID = 6
    DATA_LINE_ID = 7
    
    Z_SCORE_WARN = 3.0
    Z_SCORE_FAIL = 5.0
    
    @classmethod
    def _parse_data(cls, calling_report: str) -> dict:
        """
        :param calling_report: *.variant_calling_summary_metrics
        :return:
        """
        with open(calling_report, 'r') as in_file:
            lines = in_file.readlines()
        
        assert len(lines) > cls.DATA_LINE_ID
        
        header = lines[cls.HEADER_LINE_ID][:-1].split('\t')
        data = lines[cls.DATA_LINE_ID][:-1].split('\t')
        
        return dict(zip(header, cls._numify(data)))
    
    @classmethod
    def _numify(cls, data) -> list:
        new_data = []
        for value in data:
            try:
                value = float(value)
            except ValueError:
                value = None
            
            new_data.append(value)
        
        return new_data
    
    @classmethod
    def _parse_stats(cls, panel_stats: str) -> dict:
        stats = {}
        with open(panel_stats, 'r') as stats_file:
            for line in stats_file:
                line = line.strip()
                if not line:
                    continue
                
                key, mean, std = line.split('\t')
                stats[key] = {'mean': float(mean), 'std': float(std)}
        
        return stats
    
    @classmethod
    def _zscore2status(cls, value):
        status = Status.PASS
        if abs(value) >= cls.Z_SCORE_FAIL:
            status = Status.FAIL
        elif abs(value) >= cls.Z_SCORE_WARN:
            status = Status.WARN
        
        return status
    
    @classmethod
    def _parse_config(cls, format_filename: str) -> list:
        result = []
        with open(format_filename, 'r') as in_file:
            for line in in_file:
                line = line.strip()
                if not line:
                    continue
                
                key, name, value_type = line.split('\t')
                result.append((key, name, value_type))
        
        return result
    
    @classmethod
    def _format(cls, value_type: str):
        if value_type == 'int':
            result = '%d'
        elif value_type == 'float':
            result = '%.3f'
        elif value_type == 'pct':
            result = '%.2f%%'
        else:
            result = '%s'
        return result
    
    @classmethod
    def _value(cls, config_format: dict, value: float) -> float:
        result = value
        if config_format == 'int':
            result = int(round(result))
        if config_format == 'pct':
            result = value * 100
        return result
    
    @classmethod
    def _normal_validation(cls, name: str, value: float, value_type: str, mean: float, std: float) -> dict:
        """
        Normal is good.
        :param name:
        :param value:
        :param value_type:
        :param mean:
        :param std:
        :return:
        """
        if value is None:
            z_score = ''
            z_score_str = ''
            status = Status.NONE
            value = ''
            value_str = ''
        else:
            z_score = 0 if std == 0 else (value - mean) / std
            z_score_str = '%.2f' % z_score
            status = cls._zscore2status(z_score)
            value = value * 100 if value_type == 'pct' else value
            value_str = cls._format(value_type) % value
        
        return {
            'name': name,
            'status': status,
            'value': value,
            'value_str': value_str,
            'z_score': z_score,
            'z_score_str': z_score_str
        }
    
    @classmethod
    def _zero_validation(cls, name: str, value: float, value_type: str) -> dict:
        """
        Zero is good.
        :param name:
        :param value:
        :param value_type:
        :return:
        """
        if value is None:
            status = Status.NONE
            value = ''
            value_str = ''
        else:
            if value == 0:
                status = Status.PASS
            elif 1 <= value <= 5:
                status = Status.WARN
            else:
                status = Status.FAIL

            value = value * 100 if value_type == 'pct' else value
            value_str = cls._format(value_type) % value
        
        return {
            'name': name,
            'status': status,
            'value': value,
            'value_str': value_str,
            'z_score': '',
            'z_score_str': ''
        }
    
    @classmethod
    def parse_summary(cls, data_filename: str, config_filename: str, stats_filename: str) -> list:
        """
        https://broadinstitute.github.io/picard/picard-metric-definitions.html#CollectVariantCallingMetrics.VariantCallingDetailMetrics
        :param data_filename:
        :param config_filename:
        :param stats_filename:
        :return:
        """
        data = cls._parse_data(data_filename)  # type: dict
        stats = cls._parse_stats(stats_filename)  # type: dict
        config = cls._parse_config(config_filename)  # type: list
        
        summary = []
        for key, name, value_type in config:
            assert key in data
            value = data[key]
            if key in stats:
                mean = stats[key]['mean']
                std = stats[key]['std']
                summary.append(cls._normal_validation(name, value, value_type, mean, std))
            else:
                summary.append(cls._zero_validation(name, value, value_type))
        
        return summary


class QcReporter:
    def __init__(self, src_dir: str):
        self._src_dir = src_dir
        template_dir = self._src_dir + '/rules/shared/variant/report/summary/templates'
        env = Environment(loader=FileSystemLoader(template_dir), autoescape=True)
        self._xml_template = env.get_template('quality_report.xml')
        self._html_template = env.get_template('quality_report.html')
    
    def run(
            self,
            sample_name: str,
            panel_name: str,
            original_fastqc_r1: str,
            original_fastqc_r2: str,
            trimmed_fastqc_r1: str,
            trimmed_fastqc_r2: str,
            qualimap_txt: str,
            qualimap_html: str,
            out_xml: str,
            out_html: str,
            out_pdf: str,
            calling_report: str = None
    ):
        original_summary = FastQCParser.parse_summary_pair(original_fastqc_r1, original_fastqc_r2)
        trimmed_summary = FastQCParser.parse_summary_pair(trimmed_fastqc_r1, trimmed_fastqc_r2)
        qualimap_summary = QualimapParser.parse_summmary(qualimap_txt, qualimap_html)
        
        format_filename = '%s/rules/variant/report/calling/stats/name_format' % self._src_dir
        stats_filename = '%s/rules/variant/report/calling/stats/%s_mean_std' % (self._src_dir, panel_name)
        
        if calling_report is not None and os.path.isfile(calling_report) and os.path.isfile(stats_filename):
            calling_summary = GatkCallingParser.parse_summary(calling_report, format_filename, stats_filename)
        else:
            calling_summary = []
        
        # generate XML
        xml_string = self._xml_template.render({
            'sample_name': sample_name,
            'raw_fastqc': original_summary,
            'trimmed_fastqc': trimmed_summary,
            'qualimap': qualimap_summary,
            'calling': calling_summary
        })
        with open(out_xml, 'w') as out_file:
            out_file.write(xml_string)
        
        # generate HTML
        html_string = self._html_template.render({
            'sample_name': sample_name,
            'fastqc': zip(original_summary, trimmed_summary),
            'qualimap': qualimap_summary,
            'calling': calling_summary
        })
        with open(out_html, 'w') as out_file:
            out_file.write(html_string)
        
        # generate PDF
        html = HTML(string=html_string, encoding='utf-8')
        html.write_pdf(out_pdf)
