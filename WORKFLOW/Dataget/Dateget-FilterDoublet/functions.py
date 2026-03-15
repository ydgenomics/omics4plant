class AnalysisReporter:
    def __init__(self):
        self.report = []
    def add(self, metric, value):
        self.report.append({"Metric": metric, "Value": value})
    def save(self, filename):
        pd.DataFrame(self.report).to_csv(filename, index=False)
        print(f"报告已保存至: {filename}")

reporter = AnalysisReporter()