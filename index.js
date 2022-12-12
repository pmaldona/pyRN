const fetch = require('node-fetch');
const fs = require('fs');
const path = require('path');
const semver = require('semver');
const os = require("node:os");
const yauzl = require('yauzl');
const proc = require('child_process');


function toPromise(api) {
    return function(...args) {
        return new Promise((resolve, reject) => {
            api(...args, (error, response) => {
                if(error) {
                    return reject(error);
                }
                resolve(response);
            });
        });
    }
}

async function decompress(file) {
    let yauzlPromise = toPromise(yauzl.open);
    let zipfile = await yauzlPromise(file, {lazyEntries: true});
    let openReadStream = toPromise(zipfile.openReadStream.bind(zipfile));
    zipfile.readEntry();
    zipfile.on("entry", async (entry) => {
        if(/\/$/.test(entry.fileName)) {
            console.log(`Folder?: ${entry.fileName}`)
            if(!fs.existsSync(path.join(__dirname, 'dist', entry.fileName))) {
                fs.mkdirSync(path.join(__dirname, 'dist', entry.fileName));
            }
            zipfile.readEntry();
        }
        else {
            const zipStream = fs.createWriteStream(path.join(__dirname, 'dist', entry.fileName));
            let stream = await openReadStream(entry);
            stream.on("end", () => {
                zipfile.readEntry();
            });
            stream.pipe(zipStream);
            zipStream.on("finish", () => {
                console.log(`Updated ${entry.fileName}`);
            });
        }
    });
}

async function runMain() {
    await decompress(path.join(__dirname, 'dist', 'main2.zip'));
    let command = path.join(__dirname, 'dist', 'main', 'main')
    let args = []
    // if(isDev) {
    //     console.log("DEV-Build");
    //     command = 'python';
    //     args = ['-u','main.py']
    // }
    //let child = proc.spawn('python', ['-u','main.py']);
    console.log(command);
    let child = proc.spawn(command, args);

    let eel_started = false;

    child.stdout.on('data', function(data) {
        console.log('stdout: ' + data);

        if (data.toString().includes('[eel]: Start')) {
            eel_started = true;
        }
    });

    child.stderr.on('data', function(data) {
        console.log('stderr: ' + data);
    }); 
}

runMain();