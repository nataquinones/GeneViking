import flask
import wtforms
from scripts import geneviking

app = flask.Flask(__name__)


class InputForm(wtforms.Form):
    acc = wtforms.TextField(validators=[wtforms.validators.InputRequired()])
    end = wtforms.IntegerField(validators=[wtforms.validators.InputRequired()])
    start = wtforms.IntegerField(validators=[wtforms.validators.InputRequired()])
    ndist = wtforms.IntegerField(validators=[wtforms.validators.InputRequired()])


@app.route('/', methods=['GET', 'POST'])
def index():
    form = InputForm(flask.request.form)
    if flask.request.method == 'POST' and form.validate():
        acc = form.acc.data
        start = form.start.data
        end = form.end.data
        ndist = form.ndist.data

        result = geneviking.gene_viking(acc, start, end, ndist, 'tmp/geneviking_result.tsv')
        result = flask.Markup(result.to_html(classes='min'))

        return flask.render_template("result.html",
                                    form=form,
                               acc=acc,
                               start=start,
                               end=end,
                               ndist=ndist,
                               result=result)
    else:
        return flask.render_template("home.html", form=form)


@app.route('/about', methods=['GET', 'POST'])
def about():
  return flask.render_template('about.html')


@app.route('/download')
def download_tsv():
    return flask.send_file('tmp/geneviking_result.tsv',
                           mimetype='text',
                           attachment_filename='result.tsv',
                           as_attachment=True)


if __name__ == '__main__':
    app.run(debug=True)
