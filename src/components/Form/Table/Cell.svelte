<script>
    export let content = undefined;
    export let style = "border: 0px solid black; margin: 5px; padding: 0px;";
</script>
<td style={style}>
    {#if content}
        {content}
    {:else}
        <slot></slot>
    {/if}
</td>