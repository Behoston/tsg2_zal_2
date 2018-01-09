import re
import uuid

from sanic import Sanic
from sanic import response
from sanic_jinja2 import SanicJinja2

from web.model import Job

app = Sanic()
jinja = SanicJinja2(app)


@app.route('/')
async def index(request):
    return jinja.render('index.html', request, algorithms=['SSS', 'BBB'])


@app.route('/schedule', methods={'POST'})
async def schedule(request):
    job_id = save_job(
        fasta=request.files['fasta_file'][0].body.decode(),
        algorithm=request.form['algorithm'][0],
    )
    return jinja.render('scheduled.html', request, id=job_id)


def save_job(fasta, algorithm):
    fasta = serialize_fasta(fasta)
    job = Job(fasta=fasta, algorithm=algorithm)
    job.save()
    return job.id


def serialize_fasta(fasta):
    reads = []
    for line in fasta.split('\n'):
        if line:
            if line[0] == '>':
                reads.append('')
            else:
                line = re.sub(r'\s', '', line)
                reads[-1] = f'{reads[-1]}{line}'
    return '\n'.join(reads)


@app.route('/status/<id>')
async def get_results(request, id):
    id = uuid.UUID(id)
    return jinja.render('status.html', request, done=True, doing=False, id=id)


@app.route('/download/<id>')
async def download(request, id):
    id = uuid.UUID(id)

    async def streaming_fn(response):
        response.write(f'>{id}\n')

    return response.stream(streaming_fn, headers={
        'Content-Disposition': f'attachment; filename="{id}.fasta"',
    })


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8000)
