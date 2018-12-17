from layers import ExpressionLayer


class SRAExpressionLayer(ExpressionLayer):

    _metadata = ['run_to_disease', 'run_to_study']

    def __init__(self, *args, run_to_disease=None, run_to_study=None, **kwargs):
        self.run_to_disease = run_to_disease
        self.run_to_study = run_to_study
        super().__init__(*args, **kwargs)
