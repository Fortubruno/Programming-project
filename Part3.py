from flask import Flask, render_template, send_from_directory, request, redirect, send_file
from flask_wtf import FlaskForm
from wtforms import SubmitField, StringField
from Part1 import *
from Part2 import *


def create_obj(kind, seq, name):
    if kind == "DNA":
        return DNA(seq, name)
    elif kind == "mRNA":
        return mRNA(seq, name)
    return AAChain.create_chain(seq, name)


app = Flask(__name__)

app.secret_key = "Mnemosyne_Genomics_2024"

Mnemosyne_Dataset = DATASET()


@app.route("/")
def home():
    return render_template('home.html')


@app.route('/serve=<string:x>')
def serve_something(x):
    x = x.replace("-", "/")
    return send_from_directory('.', 'static/' + x)


@app.route('/download=<string:name>&<kind>&<string:ext>', methods=['POST'])
def download_file(name, kind, ext):
    if kind in ("Oligopeptide", "Protein"):
        seq = AAChain.find_dataset(filter_by_name=name)
    else:
        seq = Mnemosyne_Dataset.find(name)

    if ext == "fa":
        dwnld_file = Fasta.create_FASTA(name, seq.str_seq)
    else:
        dwnld_file = seq.stats.create_file()

    return send_file(dwnld_file, as_attachment=True,
                     download_name=f'{name}.{ext}', mimetype='text/plain')


@app.route('/Sequence=<string:name>&<string:kind>', methods=['GET', 'POST'])
def upload_file(name, kind):
    if name != "Home":
        if kind in ("AAChain", "Oligopeptide", "Protein"):
            obj = AAChain.find_dataset(filter_by_name=name)
        else:
            obj = Mnemosyne_Dataset.find(name)

    else:
        file = request.files['file'].read()
        kind, name, seq = Mnemosyne_Dataset.add_FASTA(file.decode('utf-8'))
        obj = create_obj(kind, seq, name)

        if kind == "AA_chain":
            Mnemosyne_Dataset.remove(name)
        else:
            Mnemosyne_Dataset.replace_FASTA(obj)
    hist = obj.stats.make_histogram()

    return render_template('Sequence.html', obj=obj, hist=hist)


@app.route('/transcription=<string:name>&<string:how>&<kind>', methods=['GET', 'POST'])
def transcribe(name, how, kind):

    if kind in ("Oligopeptide", "Protein", "AAChain"):
        seq = AAChain.find_dataset(filter_by_name=name)
    else:
        seq = Mnemosyne_Dataset.find(name)

    if how == "get":
        how = request.form.get("value")

    if how == "NT":
        new_obj = seq.transcribe()

    elif how == "RT":
        new_obj = seq.transcribe(reverse=True)

    else:
        Mnemosyne_Dataset.add_obj(seq.transcribe())
        Mnemosyne_Dataset.add_obj(seq.transcribe(reverse=True))
        return redirect(f'/my_sequences={name}&mRNA')

    if kind not in ("Oligopeptide", "Protein", "AAChain"):
        Mnemosyne_Dataset.add_obj(new_obj)

    return redirect('/Sequence=' + f"{new_obj.name}&{new_obj.type}")


@app.route("/translation=<string:name>")
def translation(name):
    if len(AAChain.find_dataset(name)) == 0:
        mrna = Mnemosyne_Dataset.find(name)
        mrna.translate()
        AA_dataset = AAChain.find_dataset(filter_by_mRNA=name)
        if len(AA_dataset) == 1:
            return redirect("/Sequence=" + f"{AA_dataset['Name'].values[0]}&AAChain")
    return redirect(f"/my_sequences={name}&AAChain")


class OrderForm(FlaskForm):
    what = StringField()
    asc = SubmitField('↓')
    desc = SubmitField('↑')


@app.route('/my_sequences=<string:name>&<string:kind>', methods=['POST', 'GET'])
def my_sequences(name, kind):
    form = OrderForm()
    dataset = pd.DataFrame()
    AA_dataset = pd.DataFrame()

    if name == "each":
        dataset = Mnemosyne_Dataset.find(each=True)
        AA_dataset = AAChain.find_dataset()
    elif kind == "mRNA":
        dataset = Mnemosyne_Dataset.find(name, exactly=False)
    else:
        AA_dataset = AAChain.find_dataset(filter_by_mRNA=name)

    if form.validate_on_submit():
        if name == "each":
            AA_dataset = AAChain.order_dataset(form.what.data, form.asc.data)
        else:
            if form.what.data == "Name":
                AA_dataset = AA_dataset.sort_values("Seq", ascending=form.asc.data)
            else:
                AA_dataset = AA_dataset.sort_values(form.what.data, ascending=form.asc.data)

    return render_template('my_sequences.html', form=form, dataset=dataset, AA_db=AA_dataset)


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8000, debug=True)
