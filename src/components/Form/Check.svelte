<script>
    export let labelText = "";
    export let labelStyle = "";
    export let bind;
</script>

<label style={labelStyle}>
    {labelText}
    <input type="checkbox" bind:checked={bind} style="opacity: 1; position: relative;">
</label>