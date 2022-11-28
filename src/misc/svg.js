C2S.prototype.circle = CanvasRenderingContext2D.prototype.circle;
C2S.prototype.square = CanvasRenderingContext2D.prototype.square;
C2S.prototype.triangle = CanvasRenderingContext2D.prototype.triangle;
C2S.prototype.triangleDown = CanvasRenderingContext2D.prototype.triangleDown;
C2S.prototype.star = CanvasRenderingContext2D.prototype.star;
C2S.prototype.diamond = CanvasRenderingContext2D.prototype.diamond;
C2S.prototype.roundRect = CanvasRenderingContext2D.prototype.roundRect;
C2S.prototype.ellipse_vis = CanvasRenderingContext2D.prototype.ellipse_vis;
C2S.prototype.database = CanvasRenderingContext2D.prototype.database;
C2S.prototype.arrowEndpoint = CanvasRenderingContext2D.prototype.arrowEndpoint;
C2S.prototype.circleEndpoint = CanvasRenderingContext2D.prototype.circleEndpoint;
C2S.prototype.dashedLine = CanvasRenderingContext2D.prototype.dashedLine;

export default function exportSvg(network) {
    var networkContainer = network.get_network().body.container;
    var ctx = new C2S({width: networkContainer.clientWidth, height: networkContainer.clientWidth, embedImages: true});

    var canvasProto = network.get_network().canvas.__proto__;
    var currentGetContext = canvasProto.getContext;
    canvasProto.getContext = function()
    {
        return ctx;
    }
    var svgOptions = {
        nodes: {
            shapeProperties: {
                interpolation: false //so images are not scaled svg will get full image
            },
            scaling: { label: { drawThreshold : 0} },
            font:{color:'#000000'}
        },
        edges: {
            scaling: { label: { drawThreshold : 0} }
        }
    };
    network.get_network().setOptions(svgOptions);
    network.get_network().redraw();
    network.get_network().setOptions(options);
    canvasProto.getContext = currentGetContext;
    ctx.waitForComplete(function()
        {
            var svg = ctx.getSerializedSvg();
            showSvg(svg);
        });
}

function showSvg(svg) {
    var svgBlob = new Blob([svg], {type: 'image/svg+xml'});
    openBlob(svgBlob, "network.svg");
}

function openBlob(blob, fileName) {
    if(window.navigator && window.navigator.msSaveOrOpenBlob)
    {

        //blobToDataURL(blob, function(dataurl){window.open(dataurl);});
        window.navigator.msSaveOrOpenBlob(blob,fileName);
    }
    else
    {
        var a = document.getElementById("blobLink");
        if(!a)
        {
            a = document.createElement("a");
            document.body.appendChild(a);
            a.setAttribute("id", "blobLink");
            a.style = "display: none";
        }
        var data = window.URL.createObjectURL(blob);
        a.href = data;
        a.download = fileName;
        a.click();
        setTimeout(function()
            {
                // For Firefox it is necessary to delay revoking the ObjectURL
                window.URL.revokeObjectURL(data);
            }
            , 100);
    }
}