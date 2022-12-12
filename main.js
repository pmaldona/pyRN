const { app, BrowserWindow, dialog, ipcMain } = require('electron');
const proc = require('child_process');
const path = require('path');
const isDev = require('electron-is-dev');
const fetch = require('node-fetch');
const fs = require('fs');
const yauzl = require('yauzl');
const os = require("node:os");
const semver = require('semver');
const Progress = require('node-fetch-progress');

let child = null;
let win;

function readLocalFile(path) {
    const fileString = fs.readFileSync(path, 'utf8');
    return fileString;
}

function constructTarget() {
    const platform = os.platform();

    if(platform == 'darwin') {
        const arch = os.arch();
        if(arch == 'arm64') {
            return 'darwin-arm64';
        }
        return 'darwin-intel';
    }
    else if(platform == 'linux') {
        return 'linux-x64';
    }
    else {
        return 'windows-x64';
    }
}

function getBackendName() {
    const target = constructTarget();
    return `python_backend-${target}`;
}

function getAsset(name, assets) {
    let asset = null;
    assets.forEach(iterationAsset => {
        if(iterationAsset["name"].includes(name)) {
            asset = iterationAsset;
            return;
        }
    });
    return asset;
}

function getBrowserURL(asset) {
    return asset["browser_download_url"];
}

function getDestination(filename) {
    // if(isDev) {
    //     return path.join(__dirname, 'static', filename);
    // } else {
    //     return path.join(__dirname, 'dist', 'main', 'static', filename);
    // }
    return path.join(__dirname, 'dist', 'main', 'static', filename);
}

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

async function decompress(file, destination) {
    let yauzlPromise = toPromise(yauzl.open);
    let zipfile = await yauzlPromise(file, {lazyEntries: true});
    let openReadStream = toPromise(zipfile.openReadStream.bind(zipfile));
    zipfile.readEntry();
    zipfile.on("entry", async (entry) => {
        
        if(/\/$/.test(entry.fileName)) {
            // console.log(`Folder?: ${entry.fileName}`)
            if(!fs.existsSync(path.join(destination, entry.fileName))) {
                fs.mkdirSync(path.join(destination, entry.fileName));
            }
            zipfile.readEntry();
        }
        else {
            if(!file.includes('static') && entry.fileName.includes('static')) {
                zipfile.readEntry();
            } else {
                if(fs.existsSync(path.join(destination, entry.fileName))) {
                    fs.unlinkSync(path.join(destination, entry.fileName));
                }
                const zipStream = fs.createWriteStream(path.join(destination, entry.fileName));
                let stream = await openReadStream(entry);
                stream.on("end", () => {
                    zipfile.readEntry();
                });
                stream.pipe(zipStream);
                zipStream.on("finish", () => {
                    // console.log(`Updated ${entry.fileName}`);
                });
            }
        }
    });
    zipfile.on("end", async () => {
        console.log("END OF ZIP");
        if(!file.includes('static')) {
            fs.chmodSync(path.join(__dirname, 'dist', 'main', 'main'), '755');
            await checkForUpdate();
        } else {
            startBackend();
        }
    });
}

async function checkForUpdate (){
    const version = readLocalFile(path.join(__dirname, 'VERSION'));
    win.webContents.send('update_doing', 'Checking JavaScript Verison...');
    win.webContents.send('update_progress','');
    try {
        const latestReleaseResponse = await fetch('https://api.github.com/repos/pmaldona/pyRN/releases/latest');
        const data = await latestReleaseResponse.json();

        console.log(`is release: ${data['tag_name']} grater than downloaded: ${version}?`);

        if(semver.gt(semver.clean(data['tag_name']), semver.clean(version))) {
            const asset = getAsset("bundle", data["assets"]);
            const browserURL = getBrowserURL(asset);
            const downloadURLResponse = await fetch(browserURL);
            win.webContents.send('update_doing', 'Downloading JavaScript...');
            const progress = new Progress(downloadURLResponse, { throttle: 100 })
            progress.on('progress', (p) => {
                win.webContents.send(
                    'update_progress',
                    `${Math.floor(p.progress * 100)}% - ${p.doneh}/${p.totalh} - ${p.rateh} - ${p.etah}`);
            })
            const destination = getDestination('bundle.zip');
            const downLoadStream = fs.createWriteStream(destination);

            await new Promise((resolve, reject) => {
                downloadURLResponse.body.pipe(downLoadStream);
                downloadURLResponse.body.on("error", reject);
                downLoadStream.on("finish", resolve);
            });

            await decompress(destination, getDestination(''));
            fs.writeFile(path.join(__dirname, 'VERSION'), data['tag_name'], 'utf8', function (err) {
                if (err) {
                    return console.log(err);
                }
                console.log(`UPDATED TO VERSION: ${data['tag_name']}`);
            });
        } else {
            console.log("Local version is the same or newer!");
            startBackend();
        }
    } catch {
        startBackend();
    }

    
    // if(version === data['tag_name']) {
    //     console.log('VERSIONS ARE IDENTICAL');
    // } else {
    //     fs.writeFile(path.join(__dirname, 'VERSION'), data['tag_name'], 'utf8', function (err) {
    //         if (err) {
    //             return console.log(err);
    //         }
    //         console.log(`UPDATED TO VERSION: ${data['tag_name']}`);
    //     });
    // }
}

async function downloadPythonBackend(assets) {
    const asset = getAsset(getBackendName(), assets);
    const downloadURL = getBrowserURL(asset);
    console.log(downloadURL);
    const downloadURLResponse = await fetch(downloadURL);
    const progress = new Progress(downloadURLResponse, { throttle: 100 })
    win.webContents.send('update_doing', 'Donwloading Python Backend...');
    progress.on('progress', (p) => {
        win.webContents.send(
            'update_progress',
            `${Math.floor(p.progress * 100)}% - ${p.doneh}/${p.totalh} - ${p.rateh} - ${p.etah}`);
    })
    const destination = path.join(__dirname, 'dist', 'main.zip');
    const downLoadStream = fs.createWriteStream(destination);

    await new Promise((resolve, reject) => {
        downloadURLResponse.body.pipe(downLoadStream);
        downloadURLResponse.body.on("error", reject);
        downLoadStream.on("finish", resolve);
    });

    await decompress(destination, path.join(__dirname, 'dist'));
}

async function updatePython() {
    const version = semver.clean(fs.readFileSync(path.join(__dirname, 'dist', 'main', 'pyversion', 'pyversion'), 'utf8'));
    win.webContents.send('update_doing', 'Checking Python Verison...');
    try {
        const latestReleaseResponse = await fetch('https://api.github.com/repos/pmaldona/pyRN/releases/latest');
        const data = await latestReleaseResponse.json();
        // console.log(data);
        const asset = getAsset(getBackendName(), data.assets);
        const onlineVersion = semver.clean(asset.name.split('-v')[1].split('.')[0].replaceAll('_', '.'));

        if(semver.gt(onlineVersion, version)) {
            await downloadPythonBackend(data.assets);
        } else {
            await checkForUpdate();
        }
    } catch {
        await checkForUpdate();
    }
}

async function startBackend() {
    win.webContents.send('update_progress', '');
    win.webContents.send('update_doing', 'Starting Backend...');
    let command = path.join(__dirname, 'dist', 'main', 'main')
    let args = []
    if(isDev) {
        console.log("DEV-Build");
        command = 'python';
        args = ['-u','main.py']
    }
    //let child = proc.spawn('python', ['-u','main.py']);
    console.log(`...trying starting backend with command: ${command}`);
    let child = proc.spawn(command, args);
    let eel_started = false;

    child.stdout.on('data', function(data) {
        console.log('stdout: ' + data);

        if (data.toString().includes('[eel]: Start')) {
            eel_started = true;
            sleep(500).then(() => {
                win.loadURL('http://localhost:8000/');
            });
        }
    });

    child.stderr.on('data', function(data) {
        console.log('stderr: ' + data);
    });
}

async function createWindow () {
    win = new BrowserWindow({
        width: 1024,
        height: 960,
        webPreferences: {
            enableRemoteModule: true,
            preload: path.join(app.getAppPath(), 'preload.js')
        }
    });
    win.loadFile(path.join(app.getAppPath(), 'static', 'updating.html'));
    if(!isDev) {
        await updatePython();
    } else {
        startBackend();
    }
}

// icmp.send('http://localhost:8000/', "Hey, I'm sending a message!")
//     .then(obj => {
//         console.log(obj.open ? 'Done' : 'Failed')
//     })
//     .catch(err => console.log(err));

function sleep(ms) {
    return new Promise((resolve) => {
      setTimeout(resolve, ms);
    });
}

app.whenReady().then(() => {
    createWindow()

    app.on('activate', () => {
        if (BrowserWindow.getAllWindows().length === 0) {
            createWindow()
        }
    })
})

app.on('window-all-closed', () => {
    if (child != null) {
        child.kill('SIGINT')
    }
    app.quit()
})


ipcMain.handle('OPEN_FILE', async (event) => {
    const filePath = await dialog.showOpenDialog({filters: [
        { name: 'Network', extensions: ['txt', 'xml', 'sbml'] }
    ]});
    return filePath;
});

ipcMain.handle('SAVE_FILE', async (event) => {
    const filePath = await dialog.showSaveDialog();
    return filePath;
});